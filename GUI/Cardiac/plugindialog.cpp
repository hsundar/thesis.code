#include <QtGui>

#include "interfaces.h"
#include "plugindialog.h"

PluginDialog::PluginDialog(const QString &path, const QStringList &fileNames,
                           QWidget *parent)
    : QDialog(parent)
{
    label = new QLabel;
    label->setWordWrap(true);

    QStringList headerLabels;
    headerLabels << tr("Components");

    treeWidget = new QTreeWidget;
    treeWidget->setAlternatingRowColors(false);
    treeWidget->setSelectionMode(QAbstractItemView::NoSelection);
    treeWidget->setHeaderLabels(headerLabels);
    treeWidget->header()->hide();

    okButton = new QPushButton(tr("OK"));
    okButton->setDefault(true);

    connect(okButton, SIGNAL(clicked()), this, SLOT(close()));

    QGridLayout *mainLayout = new QGridLayout;
    mainLayout->setColumnStretch(0, 1);
    mainLayout->setColumnStretch(2, 1);
    mainLayout->addWidget(label, 0, 0, 1, 3);
    mainLayout->addWidget(treeWidget, 1, 0, 1, 3);
    mainLayout->addWidget(okButton, 2, 1);
    setLayout(mainLayout);

    interfaceIcon.addPixmap(style()->standardPixmap(QStyle::SP_DirOpenIcon),
                            QIcon::Normal, QIcon::On);
    interfaceIcon.addPixmap(style()->standardPixmap(QStyle::SP_DirClosedIcon),
                            QIcon::Normal, QIcon::Off);
    featureIcon.addPixmap(style()->standardPixmap(QStyle::SP_FileIcon));

    setWindowTitle(tr("Plugin Information"));
    populateTreeWidget(path, fileNames);
}

void PluginDialog::populateTreeWidget(const QString &path,
                                      const QStringList &fileNames)
{
    if (fileNames.isEmpty()) {
        label->setText(tr("Cardiac Viewer couldn't find any plugins in the %1 "
                          "directory.")
                       .arg(QDir::convertSeparators(path)));
        treeWidget->hide();
    } else {
        label->setText(tr("Cardiac Viewer found the following plugins in the %1 "
                          "directory:")
                       .arg(QDir::convertSeparators(path)));

        QDir dir(path);

        foreach (QString fileName, fileNames) {
            QPluginLoader loader(dir.absoluteFilePath(fileName));
            QObject *plugin = loader.instance();

            QTreeWidgetItem *pluginItem = new QTreeWidgetItem(treeWidget);
            pluginItem->setText(0, fileName);
            treeWidget->setItemExpanded(pluginItem, true);

            QFont boldFont = pluginItem->font(0);
            boldFont.setBold(true);
            pluginItem->setFont(0, boldFont);

            if (plugin) {
                RegInterface *iReg = qobject_cast<RegInterface *>(plugin);
                if (iReg)
                    addItems(pluginItem, "RegInterface", iReg->algos());

                FilterInterface *iFilter =
                        qobject_cast<FilterInterface *>(plugin);
                if (iFilter)
                    addItems(pluginItem, "FilterInterface", iFilter->filters());
            }
        }
    }
}

void PluginDialog::addItems(QTreeWidgetItem *pluginItem,
                            const char *interfaceName,
                            const QStringList &features)
{
    QTreeWidgetItem *interfaceItem = new QTreeWidgetItem(pluginItem);
    interfaceItem->setText(0, interfaceName);
    interfaceItem->setIcon(0, interfaceIcon);

    foreach (QString feature, features) {
        if (feature.endsWith("..."))
            feature.chop(3);
        QTreeWidgetItem *featureItem = new QTreeWidgetItem(interfaceItem);
        featureItem->setText(0, feature);
        featureItem->setIcon(0, featureIcon);
    }
}
