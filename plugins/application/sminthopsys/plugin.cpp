/*
 * Sminthopsys - Dunnart Systems Biology plugin
 *
 * Copyright (C) 2011-2012  Monash University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software 
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
 * MA  02110-1301, USA.
 * 
 *
 * Author(s): Steven Kieffer  <http://skieffer.info>
*/

//! @file
//! Plugin that provides methods for manipulating biological pathways
//! in Systems Biology diagrams.

#include <QObject>
#include <QtGui>
#include <QtCore>
#include <QtAlgorithms>
#include <QFileInfo>
#include <QDomDocument>
#include <QSvgGenerator>
#include <QMap>
#include <QList>

#include <string>
#include <sstream>
#include <math.h>

#include "sbml/SBMLTypes.h"

#include "libdunnartcanvas/applicationplugininterface.h"
#include "libdunnartcanvas/canvas.h"
#include "libdunnartcanvas/canvasitem.h"
#include "libdunnartcanvas/undo.h"
#include "libdunnartcanvas/pluginshapefactory.h"
#include "libdunnartcanvas/canvasapplication.h"
#include "libdunnartcanvas/shapeplugininterface.h"
#include "libdunnartcanvas/shape.h"
#include "libdunnartcanvas/fileioplugininterface.h"
#include "libdunnartcanvas/connector.h"
#include "libdunnartcanvas/distribution.h"
#include "libdunnartcanvas/guideline.h"

#include "dsbspecies.h"
#include "dsbreaction.h"
#include "dsbcompartment.h"
#include "findBranchesDialog.h"
#include "cdSpeciesAlias.h"

using namespace dunnart;

#include "pdepn.h"
#include "pdunspecifiedepn.h"
#include "pdsimplechemepn.h"
#include "pdmacromolepn.h"
#include "pdnucleicepn.h"
#include "pdperturbingepn.h"
#include "pdsourcesink.h"
#include "pdcomplexepn.h"
#include "pdprocessnode.h"
#include "pdphenotypeprocessnode.h"
#include "sbgnGlyph.h"
#include "sbgnArc.h"

//! @brief  Plugin that provides methods for manipulating biological pathways
//! in systems biology diagrams.
//!
class SminthopsysPlugin : public QObject, public ApplicationPluginInterface,
        public ShapePluginInterface, public FileIOPluginInterface
{
    Q_OBJECT
    Q_INTERFACES (dunnart::ApplicationPluginInterface)
    Q_INTERFACES (dunnart::FileIOPluginInterface)
    Q_INTERFACES (dunnart::ShapePluginInterface)

    public:
        SminthopsysPlugin()
            : m_canvas_application(NULL)
        {
        }
        virtual void applicationMainWindowInitialised(
                CanvasApplication *canvasApplication)
        {
            m_canvas_application = canvasApplication;

            QAction *findBranchesAction = new QAction(tr("Find branches..."), this);
            connect(findBranchesAction, SIGNAL(triggered()), this,
                    SLOT(showFindBranchesDialog()));

            QMenu *dsbpeMenu = canvasApplication->mainWindow()->
                    menuBar()->addMenu("SB Pathways");
            dsbpeMenu->addAction(findBranchesAction);

        }
        virtual void applicationWillClose(CanvasApplication *canvasApplication)
        {
            Q_UNUSED (canvasApplication)
        }

        // ------------------------------------------------------------
        // SBGN Shapes

        QString shapesClassLabel(void) const
        {
            return "SBGN";
        }
        QStringList producableShapeTypes() const
        {
            QStringList shapes;
            shapes << "org.sbgn.pd.00UnspecifiedEPN";
            shapes << "org.sbgn.pd.01SimpleChemEPN";
            shapes << "org.sbgn.pd.02MacromolEPN";
            shapes << "org.sbgn.pd.03NucleicAcidEPN";
            shapes << "org.sbgn.pd.04PerturbingEPN";
            shapes << "org.sbgn.pd.05SourceOrSink";
            shapes << "org.sbgn.pd.06ComplexEPN";
            shapes << "org.sbgn.pd.ProcessNodeVertical";
            shapes << "org.sbgn.pd.ProcessNodeHorizontal";
            shapes << "org.sbgn.pd.UnknownProcessNode";
            shapes << "org.sbgn.pd.OmittedProcessNode";
            shapes << "org.sbgn.pd.AssociationProcessNode";
            shapes << "org.sbgn.pd.DissociationProcessNode";
            shapes << "org.sbgn.pd.PhenotypeProcessNode";
            return shapes;
        }
        ShapeObj *generateShapeOfType(QString shapeType)
        {
            if (shapeType == "org.sbgn.pd.00UnspecifiedEPN")
            {
                return new UnspecifiedEPN("LABEL", false);
            }
            else if (shapeType == "org.sbgn.pd.01SimpleChemEPN")
            {
                return new SimpleChemEPN("LABEL", false, false);
            }
            else if (shapeType == "org.sbgn.pd.02MacromolEPN")
            {
                return new MacromolEPN("LABEL", false, "", false);
            }
            else if (shapeType == "org.sbgn.pd.03NucleicAcidEPN")
            {
                return new NucleicAcidEPN("LABEL", false, "", false);
            }
            else if (shapeType == "org.sbgn.pd.04PerturbingEPN")
            {
                return new PerturbingEPN("LABEL", false);
            }
            else if (shapeType == "org.sbgn.pd.05SourceOrSink")
            {
                return new SourceOrSink();
            }
            else if (shapeType == "org.sbgn.pd.06ComplexEPN")
            {
                return new ComplexEPN("LABEL", false, "", false);
            }
            else if (shapeType == "org.sbgn.pd.ProcessNodeVertical")
            {
                return new ProcessNode(Qt::Vertical, ProcessNode::PROCESS);
            }
            else if (shapeType == "org.sbgn.pd.ProcessNodeHorizontal")
            {
                return new ProcessNode(Qt::Horizontal, ProcessNode::PROCESS);
            }
            else if (shapeType == "org.sbgn.pd.UnknownProcessNode")
            {
                return new ProcessNode(Qt::Vertical, ProcessNode::UNCERTAIN);
            }
            else if (shapeType == "org.sbgn.pd.OmittedProcessNode")
            {
                return new ProcessNode(Qt::Horizontal, ProcessNode::OMITTED);
            }
            else if (shapeType == "org.sbgn.pd.AssociationProcessNode")
            {
                return new ProcessNode(Qt::Horizontal, ProcessNode::ASSOCIATION);
            }
            else if (shapeType == "org.sbgn.pd.DissociationProcessNode")
            {
                return new ProcessNode(Qt::Horizontal, ProcessNode::DISSOCIATION);
            }
            else if (shapeType == "org.sbgn.pd.PhenotypeProcessNode")
            {
                return new PhenotypeProcessNode("LABEL", false);
            }
            return NULL;
        }

        // --------------------------------------------------------
        // SBML / SBGN-ML FileIO

        QStringList saveableFileExtensions(void) const
        {
            QStringList fileTypes;
            fileTypes << "sbml";
            return fileTypes;
        }
        QStringList loadableFileExtensions(void) const
        {
            QStringList fileTypes;
            fileTypes << "xml";
            fileTypes << "sbml";
            fileTypes << "sbgn";
            return fileTypes;
        }
        QString fileExtensionDescription(const QString& extension) const
        {
            if (extension == "sbml")
            {
                return "Systems Biology Markup Language (sbml)";
            }
            else if (extension == "xml")
            {
                return "Systems Biology Markup Language (xml)";
            }
            else if (extension == "sbgn")
            {
                return "Systems Biology Graphical Notation Markup Language";
            }
            return QString();
        }
        bool saveDiagramToFile(Canvas *canvas, const QFileInfo& fileInfo,
                QString& errorMessage)
        {
            Q_UNUSED(canvas);
            Q_UNUSED(fileInfo);
            Q_UNUSED(errorMessage);
            // TODO
            return false;
        }
        // TODO: If not using this method, then delete it:
        static QString nodeToString(const QDomNode& node)
        {
            QString nodeString;
            QTextStream nodeTextStream(&nodeString);
            node.save(nodeTextStream, 4);

            return nodeString;
        }

        /**
         * Returns first child of parent whose tag is as given.
         * Returns parent itself if no such child is found.
         */
        QDomNode getFirstChildOfTag(QDomNode parent, QString tag)
        {
            QDomNode first = parent;
            QDomNode node = parent.firstChild();
            while (!node.isNull()) {
                QString t = node.nodeName();
                if (t==tag) {
                    first = node; break;
                }
                node = node.nextSibling();
            }
            return first;
        }

        QList<CDSpeciesAlias*> buildSpeciesAliases(QDomNode listOfSpeciesAliases)
        {
            QList<CDSpeciesAlias*> sas;
            QDomNode sa = listOfSpeciesAliases.firstChild();
            while (!sa.isNull()) {
                CDSpeciesAlias *cdsa = new CDSpeciesAlias(sa);
                sas.append(cdsa);
                sa = sa.nextSibling();
            }
            return sas;
        }

        QList<QDomNode> getTopLevelAnnotations(QDomNode node)
        {
            QList<QDomNode> tlas;
            qDebug() << "writing all top-level tags";
            node = node.firstChild();
            while (!node.isNull()) {
                //node = node.firstChild();
                QString tag = node.nodeName();
                qDebug() << tag;
                if (tag == "annotation") {
                    tlas.append(node);
                }
                node = node.nextSibling();
            }
            qDebug() << "done writing all top-level tags";
            return tlas;
        }

        bool loadDiagramFromFile(Canvas *canvas, const QFileInfo& fileInfo,
                QString& errorMessage)
        {
            QString ext = fileInfo.suffix();
            if (ext == "sbml" || ext == "xml")
            {
                return loadDiagramFromSBML(canvas, fileInfo, errorMessage);
            }
            else if (ext == "sbgn")
            {
                return loadDiagramFromSBGN(canvas, fileInfo, errorMessage);
            }
            else
            {
                errorMessage = tr("File does not have recognized extension.");
                return false;
            }
        }

        bool loadDiagramFromSBGN(Canvas *canvas, const QFileInfo& fileInfo,
                QString& errorMessage)
        {
            QString filename = fileInfo.absoluteFilePath();
            QDomDocument dom(filename);
            QFile file(filename);

            if (!file.open(QIODevice::ReadOnly))
            {
                errorMessage = tr("File could not be opened for reading.");
                return false;
            }

            QString parsingError;
            int errorLine;
            int errorColumn;
            if (!dom.setContent(&file, true, &parsingError, &errorLine, &errorColumn))
            {
                file.close();
                errorMessage = tr("Error reading SBGN: %1:%2: %3").arg(errorLine).
                        arg(errorColumn).arg(parsingError);
                return false;
            }
            file.close();

            QDomElement sbgnTag = dom.documentElement();
            QDomNode mapTag = sbgnTag.firstChild();

            //QList<SBGNGlyph*> sbgnGlyphs;
            QMap<QString, SBGNGlyph*> sbgnGlyphs;
            QList<SBGNArc*> sbgnArcs;

            QDomNode child = mapTag.firstChild();
            while(!child.isNull() && child.nodeName() == "glyph")
            {
                SBGNGlyph *glyph = new SBGNGlyph(child);
                //qDebug() << glyph->toString() << "\n";
                sbgnGlyphs.insert(glyph->id(), glyph);
                child = child.nextSibling();
            }
            while(!child.isNull() && child.nodeName() == "arc")
            {
                SBGNArc *arc = new SBGNArc(child);
                //qDebug() << arc->toString() << "\n";
                sbgnArcs.append(arc);
                child = child.nextSibling();
            }
            // test:
            //qDebug() << "number of glyph nodes: " << sbgnGlyphs.length();
            //qDebug() << "number of arc nodes: " << sbgnArcs.length();

            // Now add shapes and connectors to the canvas.

            // Shapes
            foreach(SBGNGlyph *glyph, sbgnGlyphs.values())
            {
                ShapeObj *shape = glyph->shape();
                if (shape != NULL) {
                    canvas->addItem(shape);
                }
            }

            // Connectors
            foreach(SBGNArc *arc, sbgnArcs)
            {
                QString srcId = arc->srcId();
                // Chop off .1 or .2
                if (srcId.right(2).left(1) == ".") {
                    srcId = srcId.left(srcId.length() - 2);
                }
                QString tgtId = arc->tgtId();
                // Chop off .1 or .2
                if (tgtId.right(2).left(1) == ".") {
                    tgtId = tgtId.left(tgtId.length() - 2);
                }
                SBGNGlyph *src = sbgnGlyphs.value(srcId);
                SBGNGlyph *tgt = sbgnGlyphs.value(tgtId);
                // Tell the glyphs that they are neighbours.
                src->addNeighbour(tgt);
                tgt->addNeighbour(src);
                // Create a connector.
                Connector *conn = new Connector();
                conn->initWithConnection(src->shape(), tgt->shape());
                canvas->addItem(conn);
            }

//#define autoInferVAligns
#ifdef  autoInferVAligns
            // Vertical alignments
            // Sort by x-coord
            QMap<qreal,SBGNGlyph*> xmap;
            foreach(SBGNGlyph *g, sbgnGlyphs)
            {
                xmap.insertMulti(g->cx(),g);
            }
            // Partition into equivalence classes by x-coord
            qreal eps = 1; // tolerance
            QList< QList<SBGNGlyph*> > classes;
            QList<qreal> keys = xmap.uniqueKeys();
            int N = keys.length();
            QList<SBGNGlyph*> curr;
            qreal k = keys.at(0);
            curr.append(xmap.values(k));
            qreal lastKey = k;
            for (int i = 1; i < N; i++)
            {
                qreal k = keys.at(i);
                qDebug() << "next x value: " << k;
                if (k - lastKey >= eps) {
                    // The current key is not within tolerance of the previous one.
                    // So record the current class, and then start a new one.
                    //qDebug() << "number of classes so far: " << classes.length();
                    //qDebug() << "number of x values in current component: " << curr.length();
                    classes.append(curr);
                    //qDebug() << "number of classes after appending new one: " << classes.length();
                    curr.clear();
                }
                curr.append(xmap.values(k));
                if (i == N - 1) {
                    classes.append(curr);
                }
                lastKey = k;
            }
            // Now apply a constraint to each list.
            foreach (QList<SBGNGlyph*> list, classes)
            {
                //DEBUG
                QString s = "";
                s += "Elements of list:\n  ";
                foreach (SBGNGlyph *g, list) {
                    s += g->id() + ", ";
                }
                qDebug() << s;
                //

                // Align, if at least 2 in class
                if (list.length() < 2) { continue; }
                CanvasItemList items;
                foreach (SBGNGlyph *glyph, list) {
                    items.append(glyph->shape());
                }
                createAlignment(ALIGN_CENTER, items);

                qDebug() << "number of aligned glyphs: " << list.length();

                // Partition into connected components.
                QList< QSet<SBGNGlyph*> > ccs;
                QSet<SBGNGlyph*> S = list.toSet();
                QSet<SBGNGlyph*> R = list.toSet();
                while (!S.isEmpty()) {
                    SBGNGlyph *glyph = S.toList().first();
                    QSet<SBGNGlyph*> C;
                    glyph->getRestrConnComp(R,C);
                    ccs.append(C);
                    S = S.subtract(C);
                }
                qDebug() << "number of connected components: " << ccs.size();
                // Create a distribution constraint for CCs of 3 or more elements.
                foreach (QSet<SBGNGlyph*> cc, ccs)
                {
                    if (cc.size() < 3) { continue; }
                    // Sort by y-coord and compute average separation in y-dimension.
                    QMap<qreal,SBGNGlyph*> ymap;
                    foreach(SBGNGlyph *g, cc)
                    {
                        ymap.insert(g->cy(),g);
                    }
                    QList<qreal> yVals = ymap.keys();
                    int M = yVals.length();
                    qreal dyAvg = ( yVals.at(M-1) - yVals.at(0) ) / (M-1);
                    // Create distribution.
                    CanvasItemList items;
                    foreach (SBGNGlyph *g, cc)
                    {
                        items.append(g->shape());
                    }
                    dtype type = DIST_MIDDLE;
                    bool preserveOrder = false;
                    Distribution *dist = createDistribution(NULL, type, items, preserveOrder);
                    dist->setSeparation(dyAvg);
                }
            }
#endif

            // Turn on automatic graph layout.
            //canvas->setOptAutomaticGraphLayout(true);

            return true;
        }

        bool xcoordOrder(SBGNGlyph *g1, SBGNGlyph *g2)
        {
            return g1->cx() < g2->cx();
        }

        bool ycoordOrder(SBGNGlyph *g1, SBGNGlyph *g2)
        {
            return g1->cy() < g2->cy();
        }

        bool loadDiagramFromSBML(Canvas *canvas, const QFileInfo& fileInfo,
                QString& errorMessage)
        {
            QString filename = fileInfo.absoluteFilePath();
            SBMLDocument *doc = readSBML(filename.toStdString().c_str());

            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            // Let's try to parse it as generic xml, so we can get the
            // CellDesigner annotation tags, if they are present.
            QDomDocument dom(filename);
            QFile file(filename);

            if (!file.open(QIODevice::ReadOnly))
            {
                errorMessage = tr("File could not be opened for reading.");
                return false;
            }

            QString parsingError;
            int errorLine;
            int errorColumn;
            if (!dom.setContent(&file, true, &parsingError, &errorLine, &errorColumn))
            {
                file.close();
                errorMessage = tr("Error reading SBML: %1:%2: %3").arg(errorLine).
                        arg(errorColumn).arg(parsingError);
                return false;
            }
            file.close();

            QDomElement sbmlTag = dom.documentElement();
            QDomNode modelTag = sbmlTag.firstChild();
            //QDomNode cdAnno = getFirstChildOfTag(modelTag, "annotation");
            QDomNode cdAnno = modelTag.namedItem("annotation");
            //QList<QDomNode> tlas = getTopLevelAnnotations(modelTag);

            //QDomNode losa = getFirstChildOfTag(cdAnno, "celldesigner:listOfSpeciesAliases");
            QDomNode losa = cdAnno.namedItem("listOfSpeciesAliases");

            //qDebug() << losa.nodeName();
            QList<CDSpeciesAlias*> cdsas = buildSpeciesAliases(losa);

            //qDebug() << "Found " << tlas.length() << " top-level annotations.";
            //getTopLevelAnnotations(tlas.at(0));


            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            // Check for errors.
            unsigned int errors = doc->getNumErrors();
            if (errors > 0)
            {
                std::stringstream ss;
                ss << "Error reading SBML.\n";
                doc->printErrors(ss);
                errorMessage = QString(ss.str().c_str());
                return false;
            }

            // Get the species and reactions.
            Model *model = doc->getModel();
            ListOfSpecies *los = model->getListOfSpecies();
            unsigned int numSpecies = los->size();
            Species *spec;
            std::string id;

            // Build a map from species id's to internal objects representing those species.
            QMap<QString, DSBSpecies*> speciesMap;
            // Also a map from compartment names to DSBCompartment objects.
            QMap<QString, DSBCompartment*> compMap;
            // Some species and/or reactions (usually reactions) may not have an
            // explicitly declared compartment. Prepare lists for them.
            QList<DSBSpecies*> homelessSpecies;
            QList<DSBReaction*> homelessReacs;

            for (unsigned int i = 0; i < numSpecies; i++)
            {
                spec = los->get(i);
                id = spec->getId();
                DSBSpecies *dsbspec = new DSBSpecies(spec);
                dsbspec->setCanvas(canvas);

                // Save it in the species map.
                speciesMap.insert(QString(id.c_str()), dsbspec);

                // Get compartment name.
                QString compName = dsbspec->getCompartmentName();
                // If no name...
                if (compName.size() == 0)
                {
                    // ...add to homeless list.
                    homelessSpecies.append(dsbspec);
                }
                else // it does have a name
                {
                    // If this name has not yet been seen...
                    if (!compMap.contains(compName))
                    {
                        // ...create new entry for it.
                        DSBCompartment *comp = new DSBCompartment(compName);
                        compMap.insert(compName, comp);
                        comp->setCanvas(canvas);
                    }
                    // Now add the species to its compartment.
                    DSBCompartment *comp = compMap.value(compName);
                    comp->addSpecies(dsbspec);
                    // And give it a reference to its compartment.
                    dsbspec->setCompartment(comp);
                }
            }

            // Now get the reactions.
            ListOfReactions *lor = model->getListOfReactions();
            unsigned int numReacs = lor->size();
            Reaction *reac;
            for (unsigned int i = 0; i < numReacs; i++)
            {
                reac = lor->get(i);
                id = reac->getId();
                DSBReaction *dsbreac = new DSBReaction(reac);
                dsbreac->doublyLink(speciesMap);
                dsbreac->setCanvas(canvas);

                // Get compartment name.
                QString compName = dsbreac->getCompartmentName();
                // If no name...
                if (compName.size() == 0)
                {
                    // ...add to homeless list.
                    homelessReacs.append(dsbreac);
                }
                else // it does have a name
                {
                    // If this name has not yet been seen...
                    if (!compMap.contains(compName))
                    {
                        // ...create new entry for it.
                        DSBCompartment *comp = new DSBCompartment(compName);
                        compMap.insert(compName, comp);
                        comp->setCanvas(canvas);
                    }
                    // Now add the reaction to its compartment.
                    DSBCompartment *comp = compMap.value(compName);
                    comp->addReaction(dsbreac);
                    // And give it a reference to its compartment.
                    dsbreac->setCompartment(comp);
                }
            }

            // Set all trivial clonings.
            QList<DSBCompartment*> comps = compMap.values();
            for (int i = 0; i < comps.size(); i++)
            {
                DSBCompartment *comp = comps.at(i);
                comp->setTrivialCloning();
            }
            // Turn on automatic graph layout.
            canvas->setOptAutomaticGraphLayout(true);
            // Put all compartments in a single container.
            DSBCompartment *comp = new DSBCompartment(QString("_root"));
            comp->setCanvas(canvas);
            comp->setBoundaryVisible(false);
            comp->addCompartments(compMap.values());
            // Stash homeless species and reactions in there too.
            foreach (DSBSpecies *spec, homelessSpecies)
            {
                comp->addSpecies(spec);
                spec->setCompartment(comp);
            }
            foreach (DSBReaction *reac, homelessReacs)
            {
                comp->addReaction(reac);
                reac->setCompartment(comp);
            }
            // Set trivial cloning for any homeless species.
            comp->setTrivialCloning();
            // Finally, layout and draw.
            comp->layout();
            comp->setRelPt(QPointF(0,0));
            comp->draw();
            return true;
        }

    private slots:
        void showFindBranchesDialog()
        {
            Canvas *canvas = m_canvas_application->currentCanvas();
            QMainWindow *mainWindow = m_canvas_application->mainWindow();
            FindBranchesDialog *fbDialog = new FindBranchesDialog(canvas, mainWindow);
            fbDialog->show();
            fbDialog->raise();
            fbDialog->activateWindow();
        }

    private:
        CanvasApplication *m_canvas_application;
};

Q_EXPORT_PLUGIN2(application_sminthopsys, SminthopsysPlugin)

// Because there is no header file, we need to load the MOC file here to 
// cause Qt to generate it for us.
#include "plugin.moc"

// vim: filetype=cpp ts=4 sw=4 et tw=0 wm=0 cindent
