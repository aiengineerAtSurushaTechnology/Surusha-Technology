#ifndef TEMPLATE_H
#define TEMPLATE_H

#include <QDialog>

namespace Ui {
class Template;
}

class Template : public QDialog
{
    Q_OBJECT

public:
    explicit Template(QWidget *parent = nullptr);
    ~Template();

private slots:
    void on_okButton_clicked();

    void on_cancelButton_clicked();

private:
    Ui::Template *ui;
};

#endif // TEMPLATE_H
