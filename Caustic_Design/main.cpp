#include <QtGui>
#include "window.h"


int main(int argv, char **args)
{	
    srand(1);
	QApplication app(argv, args);
	app.setApplicationName("image CCVT 2D");

    MainWindow window;



	window.show();
	return app.exec();
}
