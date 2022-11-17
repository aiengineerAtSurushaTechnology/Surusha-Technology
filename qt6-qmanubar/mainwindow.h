#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include<QFont>
#include<QPrinter>
#include<QPrintDialog>
#include<QFileDialog>
#include <QPageSetupDialog>

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
    void on_actionOpen_triggered();

    void on_actionExit_triggered();

    void on_actionPrint_triggered();

    void on_actionPage_Setup_triggered();

    void on_actionSave_As_triggered();

    void on_actionNew_Document_triggered();

    void on_actionWindow_Background_triggered();

    void on_actionCut_triggered();

    void on_actionPaste_triggered();

    void on_actionCopy_triggered();

    void on_actionUndo_triggered();

    void on_actionRedo_triggered();

    void on_actionPlain_triggered();

    void on_actionClear_triggered();

    void on_actionBold_triggered();

    void on_actionItalic_triggered();

private:
    Ui::MainWindow *ui;
    QMenu *editMenu;
    QString file_path;

};
#endif // MAINWINDOW_H
