/*
 * Sminthopsys - Dunnart Systems Biology plugin
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
 * Author(s): Steven Kieffer <http://skieffer.info>
*/

#ifndef FINDBRANCHESDIALOG_H
#define FINDBRANCHESDIALOG_H

#include <QDialog>

class QLineEdit;
class QComboBox;

namespace dunnart {

class Canvas;
class DSBClone;

class FindBranchesDialog : public QDialog
{
    Q_OBJECT

public:
    FindBranchesDialog(Canvas *canvas, QWidget *parent = 0);

private slots:
    void findBranches();
    void canvasSelectionChanged();
    void test();

private:
    void getSelectedClone();

    Canvas *m_canvas;
    QLineEdit *m_endpointEdit;
    QComboBox *m_endpointCBox;
    QComboBox *m_layoutCBox;
    QString m_endpointIDString;
    DSBClone *m_endpointClone;

};

}

#endif
