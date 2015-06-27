#include <QtGui>
#include "window.h"

int main(int argv, char **args)
{	
	QApplication app(argv, args);
	app.setApplicationName("periodic BN-OT");
	MainWindow window;

	window.show();
	return app.exec();
}
