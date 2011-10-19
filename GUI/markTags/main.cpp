#include "xform.h"

#include <QApplication>
#include <QtGui>

int main(int argc, char **argv)
{
    QApplication app(argc, argv);

    XFormWidget xformWidget(0);

    
    QStyle *arthurStyle = new ArthurStyle();
    xformWidget.setStyle(arthurStyle);

    QList<QWidget *> widgets = qFindChildren<QWidget *>(&xformWidget);
    foreach (QWidget *w, widgets)
        w->setStyle(arthurStyle);
    
    xformWidget.show();

    // QApplication::setStyle(QStyleFactory::create("Plastique"));

    return app.exec();
}

