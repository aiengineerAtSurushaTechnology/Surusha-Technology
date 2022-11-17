#include "toolbar.h"
#include "ui_toolbar.h"
#include <QLabel>
#include <QDragEnterEvent>
#include <QDragMoveEvent>
#include <QDragLeaveEvent>
#include <QDropEvent>
#include <QFileInfo>
#include <QMimeData>
#include <QMouseEvent>
#include <QApplication>
#include <QDataStream>
#include <QDrag>
#include <QPainter>

Toolbar::Toolbar(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Toolbar)
{
    setMinimumSize(250,250);
    //ui->setupUi(this);
    setAcceptDrops(true);

    QLabel * qtIcon = new QLabel(this);
    qtIcon->setPixmap(QPixmap(":/new/prefix1/Toolbar-Images/Circle.png"));
    qtIcon->move(20,10);
    qtIcon->show();
    qtIcon->setAttribute(Qt::WA_DeleteOnClose);

    QLabel * cppIcon = new QLabel(this);
    cppIcon->setPixmap(QPixmap(":/new/prefix1/Toolbar-Images/Rectangle.png"));
    cppIcon->move(200,10);
    cppIcon->show();
    cppIcon->setAttribute(Qt::WA_DeleteOnClose);

    QLabel * terminalIcon = new QLabel(this);
    terminalIcon->setPixmap(QPixmap(":/new/prefix1/Toolbar-Images/Rounded-rectangle-stroked-96.png"));
    terminalIcon->move(20,100);
    terminalIcon->show();
    terminalIcon->setAttribute(Qt::WA_DeleteOnClose);

    QLabel * sqIcon = new QLabel(this);
    sqIcon->setPixmap(QPixmap(":/new/prefix1/Toolbar-Images/Square.png"));
    sqIcon->move(210,115);
    sqIcon->show();
    sqIcon->setAttribute(Qt::WA_DeleteOnClose);

    QLabel * stIcon = new QLabel(this);
    stIcon->setPixmap(QPixmap(":/new/prefix1/Toolbar-Images/Star.png"));
    stIcon->move(20,195);
    stIcon->show();
    stIcon->setAttribute(Qt::WA_DeleteOnClose);
}

Toolbar::~Toolbar()
{
    delete ui;
}


void Toolbar::mousePressEvent(QMouseEvent *event)
{
    if(event->button() == Qt::LeftButton)
        startPos = event->pos();
    QWidget::mousePressEvent(event);
}

void Toolbar::mouseMoveEvent(QMouseEvent *event)
{
     if(event->button() & Qt::LeftButton){
         int distance = (event->pos() - startPos).manhattanLength();
         if(distance >= QApplication::startDragDistance())
         {
            QLabel * child = static_cast<QLabel *>(childAt(event->pos()));
            if(!child)
                return;
            QPixmap mPixmap = child->pixmap();

            QByteArray ba;

            QDataStream dataStream(&ba, QIODevice::WriteOnly);

            dataStream << mPixmap << QPoint(event->pos() - child->pos());

            QMimeData * mimeData = new QMimeData;
            mimeData->setData("application/x-qtcustomitem",ba);

            QDrag * drag = new QDrag(this);
            drag->setMimeData(mimeData);
            drag->setPixmap(mPixmap);


            drag->setHotSpot(event->pos() - child->pos());

            //Blur the original label
            QPixmap tempPix = mPixmap;
            QPainter painter(&tempPix);
            painter.fillRect(tempPix.rect(),QColor(127,127,127,127));

            child->setPixmap(tempPix);


            if(drag->exec(Qt::MoveAction | Qt::CopyAction,Qt::CopyAction) == Qt::MoveAction){
                //Move data
                child->close();
            }else{
                //Copy action
                child->setPixmap(mPixmap);
            }
         }
     }
}

void Toolbar::dragEnterEvent(QDragEnterEvent *event)
{
    //event->accept();

    if(event->mimeData()->hasFormat("application/x-qtcustomitem")){
        if (event->source() == this){
            event->setDropAction(Qt::MoveAction);
            event->accept();
            //event->ignore();
        }else{
            event->acceptProposedAction();
        }
    }else{
        event->ignore();
    }
}

void Toolbar::dragMoveEvent(QDragMoveEvent *event)
{
    //event->accept();

    if(event->mimeData()->hasFormat("application/x-qtcustomitem")){
        if (event->source() == this){
            event->setDropAction(Qt::MoveAction);
            event->accept();
            //event->ignore();
        }else{
            event->acceptProposedAction();
        }
    }else{
        event->ignore();
    }
}

void Toolbar::dragLeaveEvent(QDragLeaveEvent *event)
{
    //event->accept();
    QWidget::dragLeaveEvent(event);
}

void Toolbar::dropEvent(QDropEvent *event)
{

    if (event->mimeData()->hasFormat("application/x-qtcustomitem")){

        QByteArray ba = event->mimeData()->data("application/x-qtcustomitem");
        QDataStream dataStream(&ba, QIODevice::ReadOnly);

        QPixmap mPixmap;
        QPoint offset;

        dataStream >> mPixmap >> offset;

        QLabel * newLabel = new QLabel(this);
        newLabel->setPixmap(mPixmap);
        newLabel->move(event->position().toPoint() - offset);
        newLabel->setAttribute(Qt::WA_DeleteOnClose);

        if (event->source() == this){
            event->setDropAction(Qt::MoveAction);
            event->accept();
            //event->ignore();
        }else{
            event->acceptProposedAction();
        }

    }else{
        event->ignore();
    }


/*    if(event->mimeData()->hasUrls()){
        QList<QUrl> urls = event->mimeData()->urls();
        if(urls.count()>1)
            return;

        QFileInfo file(urls.at(0).toLocalFile());
        QPixmap mPixmap;
        if(isImage(file.absoluteFilePath()) && (mPixmap.load(file.absoluteFilePath()))){
            ui->label->setPixmap(mPixmap.scaled(ui->label->size()));
        }
    }
*/
}




void Toolbar::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    painter.drawRoundedRect(5,5,width()-10, height()-10,3,3);

    QWidget::paintEvent(event);
}
/*
bool Toolbar::isImage(QString fullpath)
{
    QFileInfo file(fullpath);
    return ((file.suffix() == "png") ||
            (file.suffix() == "PNG") ||
            (file.suffix() == "jpg") ||
            (file.suffix() == "JPG") ||
            (file.suffix() == "jpeg") ||
            (file.suffix() == "JPEG"));
}
*/
