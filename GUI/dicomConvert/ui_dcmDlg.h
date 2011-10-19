/********************************************************************************
** Form generated from reading ui file 'dcmDlg.ui'
**
** Created: Mon Nov 24 15:37:26 2008
**      by: Qt User Interface Compiler version 4.4.0
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_DCMDLG_H
#define UI_DCMDLG_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QTreeWidget>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_dcmDlg
{
public:
    QGridLayout *gridLayout;
    QVBoxLayout *vboxLayout;
    QPushButton *m_pScanBut;
    QPushButton *m_pSave3dBut;
    QPushButton *m_pSave4dBut;
    QSpacerItem *spacerItem;
    QLabel *m_lblPreview;
    QPushButton *m_pDirBrowseBut;
    QHBoxLayout *hboxLayout;
    QLabel *label;
    QComboBox *m_pDcmDir;
    QTreeWidget *m_pPreview;

    void setupUi(QDialog *dcmDlg)
    {
    if (dcmDlg->objectName().isEmpty())
        dcmDlg->setObjectName(QString::fromUtf8("dcmDlg"));
    dcmDlg->resize(694, 513);
    gridLayout = new QGridLayout(dcmDlg);
#ifndef Q_OS_MAC
    gridLayout->setSpacing(6);
#endif
#ifndef Q_OS_MAC
    gridLayout->setMargin(9);
#endif
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
    vboxLayout = new QVBoxLayout();
#ifndef Q_OS_MAC
    vboxLayout->setSpacing(6);
#endif
#ifndef Q_OS_MAC
    vboxLayout->setMargin(0);
#endif
    vboxLayout->setObjectName(QString::fromUtf8("vboxLayout"));
    m_pScanBut = new QPushButton(dcmDlg);
    m_pScanBut->setObjectName(QString::fromUtf8("m_pScanBut"));

    vboxLayout->addWidget(m_pScanBut);

    m_pSave3dBut = new QPushButton(dcmDlg);
    m_pSave3dBut->setObjectName(QString::fromUtf8("m_pSave3dBut"));

    vboxLayout->addWidget(m_pSave3dBut);

    m_pSave4dBut = new QPushButton(dcmDlg);
    m_pSave4dBut->setObjectName(QString::fromUtf8("m_pSave4dBut"));

    vboxLayout->addWidget(m_pSave4dBut);

    spacerItem = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

    vboxLayout->addItem(spacerItem);

    m_lblPreview = new QLabel(dcmDlg);
    m_lblPreview->setObjectName(QString::fromUtf8("m_lblPreview"));
    QSizePolicy sizePolicy(static_cast<QSizePolicy::Policy>(0), static_cast<QSizePolicy::Policy>(0));
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(m_lblPreview->sizePolicy().hasHeightForWidth());
    m_lblPreview->setSizePolicy(sizePolicy);
    m_lblPreview->setMinimumSize(QSize(128, 128));
    m_lblPreview->setMaximumSize(QSize(128, 128));
    m_lblPreview->setFrameShape(QFrame::StyledPanel);
    m_lblPreview->setScaledContents(true);
    m_lblPreview->setAlignment(Qt::AlignCenter);

    vboxLayout->addWidget(m_lblPreview);


    gridLayout->addLayout(vboxLayout, 1, 1, 1, 1);

    m_pDirBrowseBut = new QPushButton(dcmDlg);
    m_pDirBrowseBut->setObjectName(QString::fromUtf8("m_pDirBrowseBut"));

    gridLayout->addWidget(m_pDirBrowseBut, 0, 1, 1, 1);

    hboxLayout = new QHBoxLayout();
#ifndef Q_OS_MAC
    hboxLayout->setSpacing(6);
#endif
    hboxLayout->setMargin(0);
    hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
    label = new QLabel(dcmDlg);
    label->setObjectName(QString::fromUtf8("label"));
    QSizePolicy sizePolicy1(static_cast<QSizePolicy::Policy>(5), static_cast<QSizePolicy::Policy>(5));
    sizePolicy1.setHorizontalStretch(0);
    sizePolicy1.setVerticalStretch(0);
    sizePolicy1.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
    label->setSizePolicy(sizePolicy1);

    hboxLayout->addWidget(label);

    m_pDcmDir = new QComboBox(dcmDlg);
    m_pDcmDir->setObjectName(QString::fromUtf8("m_pDcmDir"));
    QSizePolicy sizePolicy2(static_cast<QSizePolicy::Policy>(7), static_cast<QSizePolicy::Policy>(0));
    sizePolicy2.setHorizontalStretch(0);
    sizePolicy2.setVerticalStretch(0);
    sizePolicy2.setHeightForWidth(m_pDcmDir->sizePolicy().hasHeightForWidth());
    m_pDcmDir->setSizePolicy(sizePolicy2);
    m_pDcmDir->setEditable(true);

    hboxLayout->addWidget(m_pDcmDir);


    gridLayout->addLayout(hboxLayout, 0, 0, 1, 1);

    m_pPreview = new QTreeWidget(dcmDlg);
    m_pPreview->setObjectName(QString::fromUtf8("m_pPreview"));

    gridLayout->addWidget(m_pPreview, 1, 0, 1, 1);


    retranslateUi(dcmDlg);

    QMetaObject::connectSlotsByName(dcmDlg);
    } // setupUi

    void retranslateUi(QDialog *dcmDlg)
    {
    dcmDlg->setWindowTitle(QApplication::translate("dcmDlg", "Dicom Import Dialog", 0, QApplication::UnicodeUTF8));
    m_pScanBut->setText(QApplication::translate("dcmDlg", "&Scan", 0, QApplication::UnicodeUTF8));
    m_pSave3dBut->setText(QApplication::translate("dcmDlg", "Save 3D", 0, QApplication::UnicodeUTF8));
    m_pSave4dBut->setText(QApplication::translate("dcmDlg", "Save 4D", 0, QApplication::UnicodeUTF8));
    m_lblPreview->setText(QApplication::translate("dcmDlg", "Preview", 0, QApplication::UnicodeUTF8));
    m_pDirBrowseBut->setText(QApplication::translate("dcmDlg", "&Browse ...", 0, QApplication::UnicodeUTF8));
    label->setText(QApplication::translate("dcmDlg", "<html><head><meta name=\"qrichtext\" content=\"1\" /></head><body style=\" white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;\"><p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Choose <span style=\" font-weight:600;\">Directory</span>:</p></body></html>", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    m_pDcmDir->setToolTip(QApplication::translate("dcmDlg", "The Dicom Directory", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP

    Q_UNUSED(dcmDlg);
    } // retranslateUi

};

namespace Ui {
    class dcmDlg: public Ui_dcmDlg {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DCMDLG_H
