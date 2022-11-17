#include "stylesheet.h"
#include "ui_stylesheet.h"
#include <QMessageBox>
#include <QTimer>

StyleSheet::StyleSheet(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::StyleSheet)
{
    ui->setupUi(this);
}

StyleSheet::~StyleSheet()
{
    delete ui;
}

void StyleSheet::on_okButton_clicked()
{
    QMessageBox::about(this, "Message", "We need to add style sheet functions and need to enable ok button.");
}


void StyleSheet::on_cancelButton_clicked()
{
    QWidget::close();
}

