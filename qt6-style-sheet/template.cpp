#include "template.h"
#include "ui_template.h"
#include <QMessageBox>

Template::Template(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Template)
{
    ui->setupUi(this);
}

Template::~Template()
{
    delete ui;
}

void Template::on_okButton_clicked()
{
    QMessageBox::about(this, "Message", "We need to add template functions and need to enable ok button.");
}


void Template::on_cancelButton_clicked()
{
    QWidget::close();
}

