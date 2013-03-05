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

using namespace dunnart;

class JSONFileIOPlugin : public QObject, public FileIOPluginInterface
{
    Q_OBJECT
        Q_INTERFACES (dunnart::FileIOPluginInterface)

    public:
        JSONFileIOPlugin()
        {
        }
        QStringList saveableFileExtensions(void) const
        {
            QStringList fileTypes;
            fileTypes << "json";
            return fileTypes;
        }
        QStringList loadableFileExtensions(void) const
        {
            QStringList fileTypes;
            fileTypes << "json";
            return fileTypes;
        }
        QString fileExtensionDescription(const QString& extension) const
        {
            if (extension == "json")
            {
                return "JavaScript Object Notation";
            }
            return QString();
        }
        bool saveDiagramToFile(Canvas *canvas, const QFileInfo& fileInfo,
                QString& errorMessage)
        {
            QString outputFilename = fileInfo.absoluteFilePath();
            QFile jsonFile(outputFilename);
            if ( ! jsonFile.open(QIODevice::WriteOnly) )
            {
                errorMessage = tr("File could not be opened for writing.");
                return false;
            }

            QString openFile = "{\n";
            jsonFile.write(openFile.toUtf8());

            QString openNodesSection = "    \"nodes\":{\n";
            jsonFile.write(openNodesSection.toUtf8());

            canvas->setRenderingForPrinting(true);
            QList<CanvasItem *> canvas_items = canvas->items();
            QList<Connector*> conns;

            QMap<int, ShapeObj*> nodes;

            QString nodesInfo = "";
            foreach (CanvasItem *item, canvas_items)
            {
                ShapeObj *sh = dynamic_cast<ShapeObj*>(item);
                if (sh)
                {   
                    nodes.insert(sh->internalId(), sh);

                    QString info = QString("        \"%1\":{\"cx\":%2,\"cy\":%3,\"w\":%4,\"h\":%5},\n");
                    info = info.arg(sh->internalId());
                    QPointF c = sh->centrePos();
                    info = info.arg(c.x()).arg(c.y());
                    QSizeF s = sh->size();
                    info = info.arg(s.width()).arg(s.height());

                    nodesInfo += info;
                    //jsonFile.write(info.toUtf8());
                }
                Connector *conn = dynamic_cast<Connector*>(item);
                if (conn) conns.append(conn);
            }
            nodesInfo = nodesInfo.left( nodesInfo.length() - 2 ) + "\n";
            jsonFile.write(nodesInfo.toUtf8());

            QString closeNodesSection = "    },\n";
            jsonFile.write(closeNodesSection.toUtf8());

            QString openEdgesSection = "    \"edges\":{\n";
            jsonFile.write(openEdgesSection.toUtf8());

            QString edgesInfo = "";
            foreach (Connector *conn, conns)
            {
                QPair<ShapeObj *, ShapeObj *> endpts = conn->getAttachedShapes();
                int srcId = endpts.first->internalId();
                int tgtId = endpts.second->internalId();

                // Write path from centre point of src to centre point of tgt.
                QString centreToCentrePath = conn->writePath();
                // Now change first and last points to src port pt and tgt port pt.
                QStringList pts = centreToCentrePath.split(" ", QString::SkipEmptyParts);
                QString newPath = "";

                ShapeObj *src = nodes.value(srcId);
                newPath += writePortPoint(src, pts.at(1));

                for (int i = 1; i < pts.length() - 1; i++)
                {
                    newPath += " "+pts.at(i);
                }

                ShapeObj *tgt = nodes.value(tgtId);
                newPath += " "+writePortPoint( tgt, pts.at( pts.length()-2 ) );

                QString info = QString("        \"%1\":{\"src\":%2,\"tgt\":%3,\"style\":\"%4\",\"pts\":\"%5\"},\n");
                info = info.arg(conn->internalId());
                info = info.arg(srcId).arg(tgtId);
                info = info.arg("ortho");
                info = info.arg(newPath);

                edgesInfo += info;
                //jsonFile.write(info.toUtf8());
            }
            edgesInfo = edgesInfo.left( edgesInfo.length() - 2 ) + "\n";
            jsonFile.write(edgesInfo.toUtf8());

            QString closeEdgesSection = "    }\n";
            jsonFile.write(closeEdgesSection.toUtf8());

            QString closeFile = "}\n";
            jsonFile.write(closeFile.toUtf8());

            jsonFile.close();
            return true;
        }

        QString writePortPoint(ShapeObj *node, QString nextPt) {
            QStringList coords = nextPt.split(",");
            float x1 = coords.at(0).toFloat();
            float y1 = coords.at(1).toFloat();

            float x0 = node->centrePos().x();
            float y0 = node->centrePos().y();
            float w = node->size().width();
            float h = node->size().height();
            float px, py;
            if (x1 > x0) {
                // Port is on right
                px = x0 + w/2.0; py = y0;
            } else if (x1 < x0) {
                // Port is on left
                px = x0 - w/2.0; py = y0;
            } else if (y1 < y0) {
                // Port is on top
                px = x0; py = y0 - h/2.0;
            } else {
                // Port is on bottom
                px = x0; py = y0 + h/2.0;
            }
            QString portPoint = QString("%1,%2").arg(px).arg(py);
            return portPoint;
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

            // TODO: Write the method!
            errorMessage = tr("Sorry, writing JSON files has been implemented, but not reading!");
            return false;
            /*
            QString parsingError;
            int errorLine;
            int errorColumn;
            if (!doc.setContent(&file, true, &parsingError, &errorLine, &errorColumn))
            {
                file.close();
                errorMessage = tr("Error reading SVG: %1:%2: %3").arg(errorLine).
                        arg(errorColumn).arg(parsingError);
                return false;
            }
            file.close();
            canvas->setSvgRendererForFile(filename);

            QDomElement root = doc.documentElement();

            for (int pass = 0; pass < PASS_LAST; ++pass)
            {
                canvas->recursiveReadSVG(root, x_dunnartNs, pass);
            }
            */

            return true;
        }

};

Q_EXPORT_PLUGIN2(fileio_json, JSONFileIOPlugin)

// Because there is no header file, we need to load the MOC file here to
// cause Qt to generate it for us.
#include "plugin.moc"
