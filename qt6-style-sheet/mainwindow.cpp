#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QFileDialog>
#include <QTimer>
#include <QDragEnterEvent>
#include <QDragMoveEvent>
#include <QDragLeaveEvent>
#include <QDropEvent>
#include <QFileInfo>
#include <QMimeData>
#include <QMessageBox>
#include "stylesheet.h"
#include "template.h"
#include "toolbar.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setAcceptDrops(true);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_actionOpen_Style_Sheet_triggered()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Style Sheet"), "C://","All files(*.*)");
    /*
    StyleSheet stylesheet;
    stylesheet.setWindowTitle("Style Sheet");
    stylesheet.setModal(true);
    stylesheet.exec();
    */
}


void MainWindow::on_actionTemplate_triggered()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Template"), "C://","All files(*.*)");
    /*
    Template templates;
    templates.setWindowTitle("Templates");
    templates.setModal(true);
    templates.exec();
    */

}


void MainWindow::on_actionSamples_triggered()
{
    QString filename = QFileDialog::getOpenFileName(this,
                                                    tr("Open Samples"),
                                                    "C://",
                                                    "All files(*)");
}


void MainWindow::on_actionClose_triggered()
{
    statusBar()->showMessage("App will be killed in 3 seconds...");
    QTimer::singleShot(3000, this, &MainWindow::quitApp);
}


void MainWindow::on_actionSave_Style_Sheet_triggered()
{
    QString fileName = QFileDialog::getSaveFileName(this,
            tr("Save Style Sheet"), "",
            tr("Address Book (*.abk);;All Files (*)"));
}


void MainWindow::on_actionSelect_All_triggered()
{

    ui->textEdit->selectAll();
}


void MainWindow::on_actionInsert_File_triggered()
{
    QString fileContent;
    //Save the file to disk
    QString filename = QFileDialog::getOpenFileName(this,  "Open File","C://");


    if(filename.isEmpty())
        return;
    QFile file(filename);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    QTextStream in(&file);
    QString line = in.readLine();

    while(!line.isNull()){
        fileContent.append(line);
        line = in.readLine();
    }
    file.close();
    ui->textEdit->clear();
    ui->textEdit->setPlainText(fileContent);
    //ui->label->setPlainText(fileContent);
}


void MainWindow::on_actionToolbar_triggered()
{

    /*QString filename = QFileDialog::getOpenFileName(this,
                                                    tr("Open Toolbar"),
                                                    "C://",
                                                    "All files(*.*)");
                                                    */

    Toolbar toolbar;
    toolbar.setWindowTitle("Toolbar");
    toolbar.setModal(true);
    toolbar.exec();
}

void MainWindow::quitApp()
{
    QApplication::quit();
    //Widget::close();
    //MainWindow::close();
}


void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
    event->accept();
}

void MainWindow::dragMoveEvent(QDragMoveEvent *event)
{
    event->accept();
}

void MainWindow::dragLeaveEvent(QDragLeaveEvent *event)
{
    event->accept();
}

void MainWindow::dropEvent(QDropEvent *event)
{
    if(event->mimeData()->hasUrls()){
        QList<QUrl> urls = event->mimeData()->urls();
        if(urls.count()>1)
            return;

        QFileInfo file(urls.at(0).toLocalFile());
        QPixmap mPixmap;

        if(isImage(file.absoluteFilePath()) && (mPixmap.load(file.absoluteFilePath()))){
        ui->label->setPixmap(mPixmap.scaled(ui->label->size()));
        }
    }

}

bool MainWindow::isImage(QString fullpath)
{
    QFileInfo file(fullpath);
    return ((file.suffix() == "png") ||
            (file.suffix() == "PNG") ||
            (file.suffix() == "jpg") ||
            (file.suffix() == "JPG") ||
            (file.suffix() == "jpeg") ||
            (file.suffix() == "JPEG"));
}

void MainWindow::on_actionClear_triggered()
{
    //ui->textEdit->setText("");
    ui->textEdit->clear();
}


void MainWindow::on_actionClear_Image_triggered()
{
    ui->label->clear();
}


void MainWindow::on_actionUndo_Text_triggered()
{
    ui->textEdit->undo();
}


void MainWindow::on_actionRedo_Text_triggered()
{
    ui->textEdit->redo();
}


void MainWindow::on_actionCut_Text_triggered()
{
    ui->textEdit->cut();
}


void MainWindow::on_actionCopy_Text_triggered()
{
    ui->textEdit->copy();
}


void MainWindow::on_actionPaste_triggered()
{
    ui->textEdit->paste();
}


void MainWindow::on_actionAbout_Software_triggered()
{
    QMessageBox::about(this, "Message",
                       "<h1 align justify style=color:blue; > This is a prototype or initial design framework of our software"
                       "which will be created for chemical discovery and drug predictio using"
                       "various Machine Learning and Deep Learning methods. This software is build by the"
                       "use of qt6 framework. For developing ML Models we used Python language and its various "
                       "libraries.</h1>"
                       "<h2 style=color:red;>Thank You.</h2>");
}


void MainWindow::on_actionAbout_Qt_triggered()
{
    QApplication::aboutQt();
}


void MainWindow::on_actionSave_Edit_File_triggered()
{
    QString nomeFile = QFileDialog::getSaveFileName(this, tr("Save text File"), "",
                                                        tr("Text File (*.txt);;C++ File (*.cpp *.h)"));
        if (nomeFile != "") {
            QFile file(nomeFile);

            if (file.open(QIODevice::ReadWrite)) {
                QTextStream stream(&file);
                stream << ui->textEdit->toPlainText();
                file.flush();
                file.close();
            }
            else {
                QMessageBox::critical(this, tr("Error"), tr("No Such File Available"));
                return;
            }
        }
}



void MainWindow::on_actionSave_Image_triggered()
{
    QImage image;
        QString imagePath = QFileDialog::getSaveFileName( this,tr("Save Image File"),"", tr("JPEG (*.jpg *.jpeg);;PNG (*.png)" ));

            QPixmap pm =ui->label->pixmap();

                if (!pm.isNull() && !imagePath.isEmpty() )
                {
                     image =pm.toImage() ;
                     image.save(imagePath);
                }
}

