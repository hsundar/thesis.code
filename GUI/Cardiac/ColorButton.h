#ifndef COLOR_BUTTON_H
#define COLOR_BUTTON_H

#include <QPushButton>

//-----------------------------------------------------------------------

class ColorButton : public QPushButton
{
    Q_OBJECT
    QColor _color;
        
public:
	ColorButton(QWidget *parent=0);
    ColorButton(const QString &name, QWidget *parent=0);
    ColorButton(QColor &c, const QString &name, QWidget *parent = 0);
    
    void   setColor(QColor &c);
    QColor color() { return(_color); }

public slots:        
    int selectColor();

signals:
    void changed(QColor &c);
};
#endif
