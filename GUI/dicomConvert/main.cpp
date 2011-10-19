#include "dcmDialog.h"

#include <QApplication>
#include <QtGui>

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    dcmDialog mainWin;
   
	// mainWin.setWindowIcon( QIcon(":/res/images/Heart.png")) ;
    mainWin.show();
   
	// QApplication::setStyle(QStyleFactory::create("WindowsXP"));
    // QApplication::setStyle(QStyleFactory::create("Plastique"));

	return app.exec();
}
