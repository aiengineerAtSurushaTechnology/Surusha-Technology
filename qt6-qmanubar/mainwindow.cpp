#include "mainwindow.h"
#include "ui_mainwindow.h"
#include<QFile>
#include<QTextStream>
#include<QFont>
#include<QFontDialog>
#include<QFileDialog>
#include<QPrintDialog>
#include<QMessageBox>
#include <QPageSetupDialog>
#include<QColorDialog>



MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    //this->setCentralWidget(ui->textEdit);
    this->setCentralWidget(ui->textEdit);
    connect(ui->actionItalic, &QAction::triggered, this, &MainWindow::on_actionItalic_triggered);
    //connect(ui->actionItalic, &QAction::triggered, this, &MyNotepad::setFontItalic);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_actionOpen_triggered()
{
    QString file_name=QFileDialog::getOpenFileName(this,"..","c://");
    QFile file(file_name);
    file_path=file_name;
    if(!file.open(QFile::ReadOnly | QFile::Text)){
        QMessageBox::warning(this,"..","file not open");
        return;
    }
    QTextStream in(&file);
    QString text=in.readAll();
    ui->textEdit->setText(text);
    file.close();
    //QMessageBox::information(this,"..","file_name");
}


void MainWindow::on_actionExit_triggered()
{
    QApplication::quit();
}

void MainWindow::on_actionPrint_triggered()
{
    QPrinter printer;
    QPrintDialog dialog(&printer,this);
    dialog.setWindowTitle("Print Document");
    if(ui->textEdit->textCursor().hasSelection())
    {
        //QAbstractPrintDialog::addEnabledOption(QAbstractPrintDialog::PrintDialogOption,true);
      dialog.QPrintDialog::options();
   // dialog.addEnabledOption(QAbstractPrintDialog::PrintSelection);
    }
    if(dialog.exec()!=QDialog::Accepted)
    {
     return;
    }
}


void MainWindow::on_actionPage_Setup_triggered()
{
   // QPrinter *printer1;
    QPageSetupDialog q;
  // QPageSetupDialog::QPageSetupDialog(QPrinter *printer1, QWidget *parent = nullptr);
}


void MainWindow::on_actionSave_As_triggered()
{
 QString file_name= QFileDialog::getSaveFileName(this,"Open the file");
 QFile file(file_name);
 file_path=file_name;
 if(!file.open(QFile::WriteOnly | QFile::Text)){
     QMessageBox::warning(this,"..","file not open");
     return;
 }
 QTextStream out(&file);
 QString text=ui->textEdit->toPlainText();
 out<<text;
 file.flush();
 file.close();
}


void MainWindow::on_actionNew_Document_triggered()
{
    QString fileName=QFileDialog::getOpenFileName(this,tr("Open file"),QString(),
    tr("Text File(*.txt)::C++ Files(*.cpp*.h)"));
    if(fileName.isEmpty())
    {
    QFile file(fileName);
    if(file.open(QIODevice::ReadOnly)){
        QMessageBox::warning(this,"..","file not open");
        return;
    }

    QTextStream in(&file);
    ui->textEdit->setText(in.readAll());
    file.close();
    }
}


void MainWindow::on_actionWindow_Background_triggered()
{
QColor color=QColorDialog::getColor(Qt::red,this,"Choose Color");
 if(color.isValid()){
  ui->textEdit->setPalette(QPalette(color));
 }
}


void MainWindow::on_actionCut_triggered()
{
    ui->textEdit->cut();
}


void MainWindow::on_actionPaste_triggered()
{
   ui->textEdit->paste();
}


void MainWindow::on_actionCopy_triggered()
{
     ui->textEdit->copy();
}


void MainWindow::on_actionUndo_triggered()
{
   ui->textEdit->undo();
}


void MainWindow::on_actionRedo_triggered()
{
 ui->textEdit->redo();
}


void MainWindow::on_actionPlain_triggered()
{

}


void MainWindow::on_actionClear_triggered()
{
    ui->textEdit->setText("");
}


void MainWindow::on_actionBold_triggered()
{

    QTextCharFormat format;
    format.setFontWeight(QFont::Bold);

    QTextCursor cursor = ui->textEdit->textCursor();
    cursor.mergeCharFormat(format);


    QTextCharFormat fmt;
    fmt.setFontWeight(QFont::Bold);
     //mergeFormatOnWordOrSelection(fmt);



}


void MainWindow::on_actionItalic_triggered()
{
bool ok;
QFont font=QFontDialog::getFont(&ok,this);
if(ok)
{
 ui->textEdit->setFont(font);
}
else return;

}

