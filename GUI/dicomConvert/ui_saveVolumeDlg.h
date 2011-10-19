/********************************************************************************
** Form generated from reading ui file 'saveVolumeDlg.ui'
**
** Created: Mon Nov 24 15:37:24 2008
**      by: Qt User Interface Compiler version 4.4.0
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_SAVEVOLUMEDLG_H
#define UI_SAVEVOLUMEDLG_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QHBoxLayout>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QTreeView>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_saveVolumeDlg
{
public:
    QVBoxLayout *vboxLayout;
    QTreeView *m_pDirView;
    QHBoxLayout *hboxLayout;
    QSpacerItem *spacerItem;
    QPushButton *m_pCreateDirBut;
    QPushButton *m_pRemoveBut;
    QHBoxLayout *hboxLayout1;
    QLabel *label_2;
    QComboBox *m_pFormat;
    QLabel *label;
    QLineEdit *m_pPrefix;
    QHBoxLayout *hboxLayout2;
    QSpacerItem *spacerItem1;
    QPushButton *pushButton_ok;
    QPushButton *pushButton_cancel;

    void setupUi(QDialog *saveVolumeDlg)
    {
    if (saveVolumeDlg->objectName().isEmpty())
        saveVolumeDlg->setObjectName(QString::fromUtf8("saveVolumeDlg"));
    saveVolumeDlg->resize(444, 369);
    vboxLayout = new QVBoxLayout(saveVolumeDlg);
#ifndef Q_OS_MAC
    vboxLayout->setSpacing(6);
#endif
#ifndef Q_OS_MAC
    vboxLayout->setMargin(9);
#endif
    vboxLayout->setObjectName(QString::fromUtf8("vboxLayout"));
    m_pDirView = new QTreeView(saveVolumeDlg);
    m_pDirView->setObjectName(QString::fromUtf8("m_pDirView"));

    vboxLayout->addWidget(m_pDirView);

    hboxLayout = new QHBoxLayout();
#ifndef Q_OS_MAC
    hboxLayout->setSpacing(6);
#endif
    hboxLayout->setMargin(0);
    hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
    spacerItem = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

    hboxLayout->addItem(spacerItem);

    m_pCreateDirBut = new QPushButton(saveVolumeDlg);
    m_pCreateDirBut->setObjectName(QString::fromUtf8("m_pCreateDirBut"));

    hboxLayout->addWidget(m_pCreateDirBut);

    m_pRemoveBut = new QPushButton(saveVolumeDlg);
    m_pRemoveBut->setObjectName(QString::fromUtf8("m_pRemoveBut"));

    hboxLayout->addWidget(m_pRemoveBut);


    vboxLayout->addLayout(hboxLayout);

    hboxLayout1 = new QHBoxLayout();
#ifndef Q_OS_MAC
    hboxLayout1->setSpacing(6);
#endif
    hboxLayout1->setMargin(0);
    hboxLayout1->setObjectName(QString::fromUtf8("hboxLayout1"));
    label_2 = new QLabel(saveVolumeDlg);
    label_2->setObjectName(QString::fromUtf8("label_2"));

    hboxLayout1->addWidget(label_2);

    m_pFormat = new QComboBox(saveVolumeDlg);
    m_pFormat->setObjectName(QString::fromUtf8("m_pFormat"));

    hboxLayout1->addWidget(m_pFormat);

    label = new QLabel(saveVolumeDlg);
    label->setObjectName(QString::fromUtf8("label"));

    hboxLayout1->addWidget(label);

    m_pPrefix = new QLineEdit(saveVolumeDlg);
    m_pPrefix->setObjectName(QString::fromUtf8("m_pPrefix"));

    hboxLayout1->addWidget(m_pPrefix);


    vboxLayout->addLayout(hboxLayout1);

    hboxLayout2 = new QHBoxLayout();
#ifndef Q_OS_MAC
    hboxLayout2->setSpacing(6);
#endif
    hboxLayout2->setMargin(0);
    hboxLayout2->setObjectName(QString::fromUtf8("hboxLayout2"));
    spacerItem1 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

    hboxLayout2->addItem(spacerItem1);

    pushButton_ok = new QPushButton(saveVolumeDlg);
    pushButton_ok->setObjectName(QString::fromUtf8("pushButton_ok"));
    pushButton_ok->setDefault(true);

    hboxLayout2->addWidget(pushButton_ok);

    pushButton_cancel = new QPushButton(saveVolumeDlg);
    pushButton_cancel->setObjectName(QString::fromUtf8("pushButton_cancel"));

    hboxLayout2->addWidget(pushButton_cancel);


    vboxLayout->addLayout(hboxLayout2);


    retranslateUi(saveVolumeDlg);
    QObject::connect(m_pPrefix, SIGNAL(returnPressed()), pushButton_ok, SLOT(animateClick()));
    QObject::connect(pushButton_ok, SIGNAL(clicked()), saveVolumeDlg, SLOT(accept()));
    QObject::connect(pushButton_cancel, SIGNAL(clicked()), saveVolumeDlg, SLOT(reject()));

    QMetaObject::connectSlotsByName(saveVolumeDlg);
    } // setupUi

    void retranslateUi(QDialog *saveVolumeDlg)
    {
    saveVolumeDlg->setWindowTitle(QApplication::translate("saveVolumeDlg", "Save Dicom Volume", 0, QApplication::UnicodeUTF8));
    m_pCreateDirBut->setText(QApplication::translate("saveVolumeDlg", "Create Directory", 0, QApplication::UnicodeUTF8));
    m_pRemoveBut->setText(QApplication::translate("saveVolumeDlg", "Remove", 0, QApplication::UnicodeUTF8));
    label_2->setText(QApplication::translate("saveVolumeDlg", "<html><head><meta name=\"qrichtext\" content=\"1\" /></head><body style=\" white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;\"><p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">Format</span></p></body></html>", 0, QApplication::UnicodeUTF8));
    m_pFormat->insertItems(0, QStringList()
     << QApplication::translate("saveVolumeDlg", "Analyze", 0, QApplication::UnicodeUTF8)
     << QApplication::translate("saveVolumeDlg", "MetaIO", 0, QApplication::UnicodeUTF8)
    );
    label->setText(QApplication::translate("saveVolumeDlg", "<html><head><meta name=\"qrichtext\" content=\"1\" /></head><body style=\" white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;\"><p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">Prefix</span></p></body></html>", 0, QApplication::UnicodeUTF8));
    pushButton_ok->setText(QApplication::translate("saveVolumeDlg", "Ok", 0, QApplication::UnicodeUTF8));
    pushButton_cancel->setText(QApplication::translate("saveVolumeDlg", "Cancel", 0, QApplication::UnicodeUTF8));
    Q_UNUSED(saveVolumeDlg);
    } // retranslateUi

};

namespace Ui {
    class saveVolumeDlg: public Ui_saveVolumeDlg {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SAVEVOLUMEDLG_H
