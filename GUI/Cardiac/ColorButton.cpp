#include "ColorButton.h"

#include <QPalette>
#include <QColorDialog>

#define COLORHEIGHT 20

//********************************************************************************/
ColorButton::ColorButton(QWidget *parent)
    : QPushButton(parent)
{
    // setFixedHeight(COLORHEIGHT);
    // setMaximumWidth(COLORHEIGHT * 2);
    
    QColor c(Qt::black);
    setColor(c);

    connect(this, SIGNAL(clicked()), SLOT(selectColor()));
}

//********************************************************************************/

ColorButton::ColorButton(const QString &name, QWidget *parent)
    : QPushButton(name, parent)
{
    setFixedHeight(COLORHEIGHT);
    setMaximumWidth(COLORHEIGHT * 2);
    
    QColor c(Qt::black);
    setColor(c);

    connect(this, SIGNAL(clicked()), SLOT(selectColor()));
}

//********************************************************************************/

ColorButton::ColorButton(QColor &c, const QString &name, QWidget *parent)
    : QPushButton(name, parent)
{
    setFixedHeight(COLORHEIGHT);
    setMaximumWidth(COLORHEIGHT * 2);
    setColor(c);
    
    connect(this, SIGNAL(clicked()), SLOT(selectColor()));
}

//********************************************************************************/

int
ColorButton::selectColor()
{
    QColor c(QColorDialog::getColor(_color));
    if(!c.isValid()) return(FALSE);

    QPalette pal(palette());
    pal.setColor(QPalette::Button, c);
    setPalette(pal);

    _color = c;
    emit changed(_color);
    return(TRUE);
}

//********************************************************************************/

void
ColorButton::setColor(QColor &c)
{
    QPalette pal(palette());
    pal.setColor(QPalette::Button, c);
    setPalette(pal);
    _color = c;
}
