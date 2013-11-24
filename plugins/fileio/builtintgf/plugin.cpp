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

            QMap<QString,ShapeObj*> shapesByNum;
            PluginShapeFactory *shapeFactory = sharedPluginShapeFactory();
            int n = 0;
            canvas->stop_graph_layout();
            while (!file.atEnd()) {
                QByteArray line = file.readLine();
                //qDebug() << line.length();
                QString s = QString(line.data()).trimmed();
                //qDebug() << s;
                QStringList parts = s.split(" ");
                //qDebug() << QString("%1 parts").arg(parts.length());
                if (parts.length() == 1) {
                    QString p = parts.at(0);
                    if (p.at(0)==QString("#").at(0)) continue;
                    ShapeObj *s = shapeFactory->createShape("org.dunnart.shapes.rect");
                    shapesByNum.insert(p,s);
                    double x = n%4<2 ? n : -n;
                    double y = n%2<1 ? n : -n;
                    s->setPosAndSize(QPointF(x,y),QSizeF(50,50));
                    canvas->addItem(s);
                    n++;
                } else if (parts.length() == 2) {
                    QString p = parts.at(0), q = parts.at(1);
                    ShapeObj *s = shapesByNum.value(p);
                    ShapeObj *t = shapesByNum.value(q);
                    Connector *c = new Connector();
                    c->initWithConnection(s,t);
                    canvas->addItem(c);
                } else {
                    // This should not happen!
                }
            }
            file.close();
            canvas->fully_restart_graph_layout();
            return true;

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

Q_EXPORT_PLUGIN2(fileio_builtintgf, BuiltinTGFFileIOPlugin)

// Because there is no header file, we need to load the MOC file here to
// cause Qt to generate it for us.
#include "plugin.moc"
