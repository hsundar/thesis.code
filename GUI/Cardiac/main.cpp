#include "mainwindow.h"
#include "arthurwidgets.h"

#include <QApplication>
#include <QtGui>

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    MainWindow mainWin;
   /*
    QStyle *arthurStyle = new ArthurStyle();
    mainWin.setStyle(arthurStyle);
    QList<QWidget *> widgets = qFindChildren<QWidget *>(&mainWin);
    foreach (QWidget *w, widgets)
        w->setStyle(arthurStyle);
   */
    mainWin.setWindowIcon( QIcon(":/res/images/Heart.png")) ;
    mainWin.show();
    QApplication::setStyle(QStyleFactory::create("WindowsXP"));
    // QApplication::setStyle(QStyleFactory::create("Cleanlooks"));
		// QApplication::setStyle(QStyleFactory::create("Plastique"));
		/*
		QFile file(":/res/cardiac.qss");
		file.open(QFile::ReadOnly);
		QString styleSheet = QLatin1String(file.readAll());

		qApp->setStyleSheet(styleSheet);
			*/

		return app.exec();
}
