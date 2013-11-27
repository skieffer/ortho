/*
 * Dunnart - Constraint-based Diagram Editor
 *
 * Copyright (C) 2011  Monash University
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
 * Author(s): Steve Kieffer  <http://skieffer.info/>
*/

#include <QPair>
#include <QtGui>
#include <QObject>
#include <QFileInfo>
#include <QDomDocument>
#include <QSvgGenerator>

#include "libdunnartcanvas/fileioplugininterface.h"
#include "libdunnartcanvas/canvas.h"
#include "libdunnartcanvas/canvasitem.h"
#include "libdunnartcanvas/shape.h"
#include "libdunnartcanvas/connector.h"

#include "libdunnartcanvas/pluginshapefactory.h"

using namespace dunnart;

class BuiltinTGFFileIOPlugin : public QObject, public FileIOPluginInterface
{
    Q_OBJECT
        Q_INTERFACES (dunnart::FileIOPluginInterface)

    public:
        BuiltinTGFFileIOPlugin()
        {
        }
        QStringList saveableFileExtensions(void) const
        {
            QStringList fileTypes;
            fileTypes << "tgf";
            return fileTypes;
        }
        QStringList loadableFileExtensions(void) const
        {
            QStringList fileTypes;
            fileTypes << "tgf";
            return fileTypes;
        }
        QString fileExtensionDescription(const QString& extension) const
        {
            if (extension == "tgf")
            {
                return "Trivial Graph Format (TGF)";
            }
            return QString();
        }
        bool saveDiagramToFile(Canvas *canvas, const QFileInfo& fileInfo,
                QString& errorMessage)
        {
            QString outputFilename = fileInfo.absoluteFilePath();
            QFile tgfFile(outputFilename);
            if ( ! tgfFile.open(QIODevice::WriteOnly) )
            {
                errorMessage = tr("File could not be opened for writing.");
                return false;
            }

            canvas->setRenderingForPrinting(true);
            QList<CanvasItem *> canvas_items = canvas->items();
            QMap<ShapeObj*,int> shapeIDs;
            QList<Connector*> conns;
            int i = 0;
            foreach (CanvasItem *item, canvas_items)
            {
                ShapeObj *sh = dynamic_cast<ShapeObj*>(item);
                if (sh)
                {
                    QString line = QString::number(i);
                    QString label = sh->getLabel();
                    if (label.length() > 0) line += " "+label;
                    line += "\n";
                    tgfFile.write(line.toUtf8());
                    shapeIDs.insert(sh,i);
                    i++;
                }
                Connector *conn = dynamic_cast<Connector*>(item);
                if (conn) conns.append(conn);
            }

            QString startEdgeSection = "#\n";
            tgfFile.write(startEdgeSection.toUtf8());

            foreach (Connector *conn, conns)
            {
                QPair<ShapeObj *, ShapeObj *> endpts = conn->getAttachedShapes();
                int srcIdx = shapeIDs.value(endpts.first);
                int tgtIdx = shapeIDs.value(endpts.second);
                QString srcStr = QString::number(srcIdx);
                QString tgtStr = QString::number(tgtIdx);
                QString edgeLine = srcStr + " " + tgtStr + "\n";
                tgfFile.write(edgeLine.toUtf8());
            }

            tgfFile.close();
            return true;
        }

        bool loadDiagramFromFile(Canvas *canvas, const QFileInfo& fileInfo,
                QString& errorMessage)
        {
            QString filename = fileInfo.absoluteFilePath();
            QDomDocument doc(filename);
            QFile file(filename);
            if (!file.open(QIODevice::ReadOnly))
            {
                errorMessage = tr("File could not be opened for reading.");
                return false;
            }
            // Read the node and edge declarations in the file.
            QTextStream in(&file);
            bool readingNodes = true;
            QList<QString> nodes;
            QList< QPair<QString,QString> > edges;
            while (!in.atEnd())
            {
                QString line = in.readLine();
                line = line.trimmed();
                if (line.length()==0) continue;
                if (line.startsWith("#")) {
                    readingNodes = false;
                } else if (readingNodes) {
                    // Record node
                    nodes.append(line);
                } else {
                    // Record edge
                    QStringList ends = line.split(" ");
                    edges.append(QPair<QString,QString>(ends.at(0), ends.at(1)));
                }
            }
            qDebug() << nodes;
            qDebug() << edges;
            // Build the nodes.
            QMap<QString, ShapeObj*> map;
            PluginShapeFactory *factory = sharedPluginShapeFactory();
            QString type("org.dunnart.shapes.rect");
            foreach (QString node, nodes) {
                ShapeObj *shape = factory->createShape(type);
                map.insert(node, shape);
                shape->setSize(QSizeF(30,30));
                shape->setCentrePos(QPointF(0,0));
                QUndoCommand *cmd = new CmdCanvasSceneAddItem(canvas, shape);
                canvas->currentUndoMacro()->addCommand(cmd);
            }
            // Build the edges.
            typedef QPair<QString,QString> Edge;
            foreach ( Edge edge, edges) {
                Connector *conn = new Connector();
                conn->initWithConnection(map.value(edge.first), map.value(edge.second));
                conn->setDirected(true);
                QUndoCommand *cmd = new CmdCanvasSceneAddItem(canvas, conn);
                canvas->currentUndoMacro()->addCommand(cmd);
            }
            return true;
        }

};

Q_EXPORT_PLUGIN2(fileio_builtintgf, BuiltinTGFFileIOPlugin)

// Because there is no header file, we need to load the MOC file here to
// cause Qt to generate it for us.
#include "plugin.moc"
