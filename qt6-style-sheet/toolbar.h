#ifndef TOOLBAR_H
#define TOOLBAR_H

#include <QDialog>

namespace Ui {
class Toolbar;
}

class Toolbar : public QDialog
{
    Q_OBJECT

public:
    explicit Toolbar(QWidget *parent = nullptr);
    ~Toolbar();

private:
    Ui::Toolbar *ui;

    // QWidget interface
protected:
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void dragEnterEvent(QDragEnterEvent *event);
    void dragMoveEvent(QDragMoveEvent *event);
    void dragLeaveEvent(QDragLeaveEvent *event);
    void dropEvent(QDropEvent *event);

private:
    QPoint startPos;
    //bool isImage(QString fullpath);

    // QWidget interface
protected:
    void paintEvent(QPaintEvent *event);
};

#endif // TOOLBAR_H
