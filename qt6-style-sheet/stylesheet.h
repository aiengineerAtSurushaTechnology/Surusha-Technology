#ifndef STYLESHEET_H
#define STYLESHEET_H

#include <QDialog>

namespace Ui {
class StyleSheet;
}

class StyleSheet : public QDialog
{
    Q_OBJECT

public:
    explicit StyleSheet(QWidget *parent = nullptr);
    ~StyleSheet();


private slots:
    void on_okButton_clicked();

    void on_cancelButton_clicked();

private:
    Ui::StyleSheet *ui;
};

#endif // STYLESHEET_H
