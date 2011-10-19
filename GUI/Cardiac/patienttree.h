#ifndef PATIENT_TREE_H
#define PATIENT_TREE_H

#include <QDomDocument>
#include <QHash>
#include <QIcon>
#include <QTreeWidget>

#include <QDir>

class patientTree : public QTreeWidget
{
    Q_OBJECT

public:
    patientTree(QWidget *parent = 0);

    void setThumbnailDirectory (QDir dir) { tnailDir = dir; };
    bool read(QIODevice *device);
    bool write(QIODevice *device);
	QDomElement initNewRoot();
	QDomDocument* getDocument() { return &domDocument; };

private slots:
    void updateDomElement(QTreeWidgetItem *item, int column);

private:
    void parsePatientElement(const QDomElement &element,
                            QTreeWidgetItem *parentItem = 0);
    QTreeWidgetItem *createItem(const QDomElement &element,
                                QTreeWidgetItem *parentItem = 0);

    QDomDocument domDocument;
    QHash<QTreeWidgetItem *, QDomElement> domElementForItem;
    
    QDir  tnailDir;
    QIcon patientIcon;
    QIcon dataIcon;
};

#endif
