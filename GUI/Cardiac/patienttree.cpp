#include <QtGui>

#include "patienttree.h"

patientTree::patientTree(QWidget *parent)
    : QTreeWidget(parent)
{
    QStringList labels;
    // labels << tr("Patient") << tr("Study") << tr("Type / Protocol");
    labels << tr("Subject") << tr("Protocol") << tr("Study ID") << tr("Series") << tr("Dimensions") << tr("Frames");

    // header()->setResizeMode(QHeaderView::Stretch);
    setHeaderLabels(labels);

    patientIcon.addPixmap(QPixmap(":/res/images/patient.png"));
    dataIcon.addPixmap ( QPixmap(":/res/images/data.png"));
}

bool patientTree::read(QIODevice *device)
{
    QString errorStr;
    int errorLine;
    int errorColumn;

    if (!domDocument.setContent(device, true, &errorStr, &errorLine,
                                &errorColumn)) {
        QMessageBox::information(window(), tr("Patient Browser"),
                                 tr("Parse error at line %1, column %2:\n%3")
                                 .arg(errorLine)
                                 .arg(errorColumn)
                                 .arg(errorStr));
        return false;
    }

    QDomElement root = domDocument.documentElement();
    if (root.tagName() != "sbia") {
        QMessageBox::information(window(), tr("Patient Browser"),
                                 tr("The file is not a SBIA patient database file."));
        return false;
    } else if (root.hasAttribute("version")
               && root.attribute("version") != "1.0") {
        QMessageBox::information(window(), tr("Patient Browser"),
                                 tr("The version of the SBIA patient database is not supported"));
        return false;
    }

    clear();

    // disconnect(this, SIGNAL(itemChanged(QTreeWidgetItem *, int)),
    //           this, SLOT(updateDomElement(QTreeWidgetItem *, int)));

    QDomElement child = root.firstChildElement("patient");
    while (!child.isNull()) {
        parsePatientElement(child);
        child = child.nextSiblingElement("patient");
    }

    // connect(this, SIGNAL(itemChanged(QTreeWidgetItem *, int)),
    //        this, SLOT(updateDomElement(QTreeWidgetItem *, int)));

    return true;
}

bool patientTree::write(QIODevice *device)
{
    const int IndentSize = 4;

    QTextStream out(device);
    domDocument.save(out, IndentSize);
    return true;
}

void patientTree::updateDomElement(QTreeWidgetItem *item, int column)
{
    QDomElement element = domElementForItem.value(item);
    
    if (!element.isNull()) {
      switch (column) {
        case 0: {
          QDomElement oldTitleElement = element.firstChildElement("name");
          QDomElement newTitleElement = domDocument.createElement("name");

          QDomText newTitleText = domDocument.createTextNode(item->text(0));
          newTitleElement.appendChild(newTitleText);

          element.replaceChild(newTitleElement, oldTitleElement);
         }
          break;
        case 1:
          if (element.tagName() == "dataset")
            element.setAttribute("protocol", item->text(1));
          break;
        case 2:
          if (element.tagName() == "dataset")
            element.setAttribute("study", item->text(2));
		  break;
		case 3:
          if (element.tagName() == "dataset")
            element.setAttribute("series", item->text(3));
      }
    }

     /*   if (column == 0) {
            QDomElement oldTitleElement = element.firstChildElement("title");
            QDomElement newTitleElement = domDocument.createElement("title");

            QDomText newTitleText = domDocument.createTextNode(item->text(0));
            newTitleElement.appendChild(newTitleText);

            element.replaceChild(newTitleElement, oldTitleElement);
        } else {
            if (element.tagName() == "bookmark")
                element.setAttribute("href", item->text(1));
        }
    } */
}

void patientTree::parsePatientElement(const QDomElement &element,
                                  QTreeWidgetItem *parentItem)
{
    QTreeWidgetItem *item = createItem(element, parentItem);

    QString title = element.firstChildElement("name").text();
    if (title.isEmpty())
        title = QObject::tr("Unknown Patient");

    item->setFlags(item->flags() | Qt::ItemIsEditable);
    item->setIcon(0, patientIcon);
    item->setText(0, title);

    bool folded = (element.attribute("folded") != "no");
    setItemExpanded(item, !folded);

    QDomElement child = element.firstChildElement();
    while (!child.isNull()) {
        if (child.tagName() == "dataset") {
            QTreeWidgetItem *childItem = createItem(child, item);

            QString proto = child.firstChildElement("protocol").text();
            if (proto.isEmpty())
                proto = QObject::tr("unknown protocol");

            QString study = child.firstChildElement("study").text();
            if (study.isEmpty())
                study = QObject::tr("?");

            QString series = child.firstChildElement("series").text();
            if (series.isEmpty())
                series = QObject::tr("?");

            QString dims = child.firstChildElement("dimensions").text();
            if (dims.isEmpty())
                dims = QObject::tr("?x?x?");

            QDomElement frames = child.firstChildElement("frames");
            QString cnt = frames.attribute("count", "1");

            // the thumbnail ...
            QString tnail = tnailDir.absoluteFilePath(frames.firstChildElement("thumbnail").text());

            childItem->setFlags(item->flags() | Qt::ItemIsEditable);
            childItem->setIcon(0, dataIcon);
            childItem->setText(1, proto);
            childItem->setText(2, study);
            childItem->setText(3, series);
            childItem->setText(4, dims);
            childItem->setText(5, cnt);
            childItem->setData(0, Qt::WhatsThisRole, QPixmap::fromImage(QImage(tnail)) );


        } else if (child.tagName() == "separator") {
            QTreeWidgetItem *childItem = createItem(child, item);
            childItem->setFlags(item->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            childItem->setText(0, QString(30, 0xB7));
        }
        child = child.nextSiblingElement();
    }
}

QTreeWidgetItem *patientTree::createItem(const QDomElement &element,
                                      QTreeWidgetItem *parentItem)
{
    QTreeWidgetItem *item;
    if (parentItem) {
        item = new QTreeWidgetItem(parentItem);
    } else {
        item = new QTreeWidgetItem(this);
    }
    domElementForItem.insert(item, element);
    return item;
}

QDomElement patientTree::initNewRoot() {
  domDocument= QDomDocument("sbia");

  QDomElement root = domDocument.createElement("sbia");
  root.setAttribute("version", "1.0");
  domDocument.appendChild(root);
  return domDocument.documentElement();
}
