#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_actionOpen_Style_Sheet_triggered();

    void on_actionTemplate_triggered();

    void on_actionSamples_triggered();

    void on_actionClose_triggered();

    void on_actionSave_Style_Sheet_triggered();

    void on_actionSelect_All_triggered();

    void on_actionInsert_File_triggered();

    void on_actionToolbar_triggered();

    void quitApp();


    void on_actionClear_triggered();

    void on_actionClear_Image_triggered();

    void on_actionUndo_Text_triggered();

    void on_actionRedo_Text_triggered();

    void on_actionCut_Text_triggered();

    void on_actionCopy_Text_triggered();

    void on_actionPaste_triggered();

    void on_actionAbout_Software_triggered();

    void on_actionAbout_Qt_triggered();

    void on_actionSave_Edit_File_triggered();

    void on_actionSave_Image_triggered();

private:
    Ui::MainWindow *ui;

protected:
    void dragEnterEvent(QDragEnterEvent *event);
    void dragMoveEvent(QDragMoveEvent *event);
    void dragLeaveEvent(QDragLeaveEvent *event);
    void dropEvent(QDropEvent *event);

private:
    bool isImage(QString fullpath);


    // QWidget interface

};
#endif // MAINWINDOW_H
