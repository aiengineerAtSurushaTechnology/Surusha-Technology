#include "dialog.h"
#include "ui_dialog.h"
#include"QGraphicsItem"

Dialog::Dialog(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::Dialog)
{
    ui->setupUi(this);
    scene =new QGraphicsScene(this);
    ui->graphicsView->setScene(scene);
    QBrush redBrush(Qt::red);
    QBrush blueBrush(Qt::blue);
    QPen blackpen(Qt::black);
blackpen.setWidth(6);

ellipse = scene->addEllipse(10,10,100,100,blackpen,redBrush);
rectangle = scene->addRect(-100,-100,50,50,blackpen,blueBrush);
rectangle->setFlag(QGraphicsItem::ItemIsMovable);
ellipse->setFlag(QGraphicsItem::ItemIsMovable);


}

Dialog::~Dialog()
{
    delete ui;
}


void Dialog::on_pushButton_clicked()
{
    ui->graphicsView->rotate(-20);
}


void Dialog::on_pushButton_2_clicked()
{
    ui->graphicsView->rotate(20);
}

