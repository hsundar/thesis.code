#include "mainwindow.h"
#include "TransMatrix.h"
#include <iostream>

MainWindow::MainWindow() {

  targetVol = NULL;
  sourceVol = NULL;
  
  sim = new oct_sim();
   
  sim->setHistBits(4);
  sim->initHists();


  QToolBar *toolBar = new QToolBar(this);
  toolBar->setFixedHeight(80);

#if QT_VERSION < 0x040000
  setDockEnabled(TornOff, true);
  setRightJustification(true);
#else
  toolBar->setAllowedAreas(Qt::TopToolBarArea | Qt::BottomToolBarArea);
#endif

  /*
   * Layout ...
   *
   * load two images .... so file browsers ...
   * then list boxes with options ...
   * finally button to test ...
   */ 

  QWidget *hBox = new QWidget(toolBar);

  // Target browser
  QLabel *label1 = new QLabel("Target", hBox);
  targetFileName = new QLineEdit(hBox);
  QToolButton *targetButton = new QToolButton(hBox);

  // Source browser
  QLabel *label2 = new QLabel("Source", hBox);
  sourceFileName = new QLineEdit(hBox);
  QToolButton *sourceButton = new QToolButton(hBox);

  // Options ...
  metricCombo = new QComboBox(hBox);
  QStringList metrices;
  metrices << "SSD" << "NCC" << "MI" << "NMI";
  metricCombo->addItems(metrices);
  metricCombo->setCurrentIndex(3);

  modesCombo = new QComboBox(hBox);
  QStringList modes;
  modes << "Trans-X" << "Trans-Y" << "Trans-Z" << "Rot-X" << "Rot-Y" << "Rot-Z" << "Rot-X,Y" << "Rot-X,Z" ;
  modesCombo->addItems(modes);

  rangeSpin = new QSpinBox(hBox);
  rangeSpin->setRange(1,50);
  rangeSpin->setValue(10);
  stepSpin = new QDoubleSpinBox(hBox);
  stepSpin->setRange(0.,5.00);
  stepSpin->setSingleStep(0.02);
  stepSpin->setValue(1.0);
  
  useOctCheck = new QCheckBox(hBox);
  useOctCheck->setChecked(true);

  runButton = new QPushButton("Run", hBox);
  runButton->setEnabled(false);

  saveButton = new QPushButton("Save", hBox);

  QHBoxLayout *layout = new QHBoxLayout(hBox);
  layout->addWidget(label1);
  layout->addWidget(targetFileName);
  layout->addWidget(targetButton); 
  layout->addWidget(label2);
  layout->addWidget(sourceFileName);
  layout->addWidget(sourceButton); 
  layout->addWidget(metricCombo);
  layout->addWidget(modesCombo);
  layout->addWidget(new QLabel("Range", hBox));
  layout->addWidget(rangeSpin);
  layout->addWidget(new QLabel("Step Size", hBox));
  layout->addWidget(stepSpin);
  layout->addWidget(useOctCheck);
  layout->addWidget(runButton);
  layout->addWidget(saveButton);
  //layout->addWidget(new QWidget(hBox), 10); // spacer);

  // Connect Signals / Slots 

  connect(targetButton, SIGNAL(clicked()), this, SLOT(onBrowseTarget()));
  connect(sourceButton, SIGNAL(clicked()), this, SLOT(onBrowseSource()));

  connect(runButton, SIGNAL(clicked()), this, SLOT(onRun()));
  connect(saveButton, SIGNAL(clicked()), this, SLOT(onSave()));


#if QT_VERSION >= 0x040000
  toolBar->addWidget(hBox);
#endif
  addToolBar(toolBar);


  m_plot = new DataPlot(this);
  setCentralWidget(m_plot);

}

void MainWindow::onBrowseTarget() {
  QString fn = QFileDialog::getOpenFileName( this, "Choose a file", "/data/hsundar/SCR", "MetaIO/Analyze Files (*.mhd *.mha *.hdr)");
  // QString fn = QFileDialog::getOpenFileName( this, "Choose a file", QDir::currentPath(), "MetaIO/Analyze Files (*.mhd *.mha *.hdr)");

  if ( !fn.isEmpty() ) {
    // delete target and load new volume ...
    if (targetVol != NULL)
      delete targetVol;
    targetVol = new Volume();

    // check format ...
    if (fn.endsWith("mhd"))
      targetVol->InitMhd(fn);
    else if ( fn.endsWith("hdr"))
      targetVol->InitAnalyze(fn);
    else {
      QMessageBox::warning(this, "Oct Similarity",
          "Cannot load Volume File.\n"
          "Please check file format.\n\n",
          "Retry", "Quit", 0, 0, 1);
      return;
    }
    // set correct paths in lineEdits ...
    targetFileName->setText(fn);
   /* if( sourceFileName->text().isEmpty() )
      sourceFileName->setText(fn);
    */

    sim->setTargetImage(targetVol);
    
    targetVol->SetMaxDepth(8);

    runButton->setEnabled(true);

    // first generate the Octree ...
    targetVol->ComputeLinearOctree();
  }
}

void MainWindow::onBrowseSource() {
  QString fn = QFileDialog::getOpenFileName( this, "Choose a file", "/data/hsundar/SCR", "MetaIO/Analyze Files (*.mhd *.mha *.hdr)");
  // QString fn = QFileDialog::getOpenFileName( this, "Choose a file", QDir::currentPath(), "MetaIO/Analyze Files (*.mhd *.mha *.hdr)");

  if ( !fn.isEmpty() ) {
    // delete target and load new volume ...
    if (sourceVol != NULL)
      delete sourceVol;
    sourceVol = new Volume();

    // check format ...
    if (fn.endsWith("mhd"))
      sourceVol->InitMhd(fn);
    else if ( fn.endsWith("hdr"))
      sourceVol->InitAnalyze(fn);
    else {
      QMessageBox::warning(this, "Oct Similarity",
          "Cannot load Volume File.\n"
          "Please check file format.\n\n",
          "Retry", "Quit", 0, 0, 1);
      return;
    }
    // set correct paths in lineEdits ...
    sourceFileName->setText(fn);
    
    sim->setSourceImage(sourceVol);

    // runButton->setEnabled(false);

  }
}

void MainWindow::onRun() {
  // will do most of what was being done in the earlier code here ...


  // Now run the loop ....

  TransMatrix mat;
  int iter = rangeSpin->value();
  double step = stepSpin->value();

  // std::cout << "range is " << iter << " and step is " << step << std::endl;
  int metric = metricCombo->currentIndex();
  int mode   = modesCombo->currentIndex(); 

  // std::cout << "Metric is " << metric << std::endl;
  // for now ... hardcode whether to use oct or not ...
  
  sim->useOctree(useOctCheck->isChecked());

  m_plot->setPlotRange(iter, step, modesCombo->currentText());
  
  // First print out the x-axis values ...
  std::cout << "xvals = [ ";
  for (int i= -iter; i<iter+1; i++)
    std::cout << i*step << ", ";
  std::cout << " ]" << std::endl << std::endl;

  std::cout << "oct_mi = [ ";

  unsigned int cnt = 0;
  double val=0;
  for (int i= -iter; i<iter+1; i++) {
    // Get the right matrix ...
    // mat = TransMatrix::RotationX(step*i);
    // mat = TransMatrix::Translation(step*i, 0., 0.);
    switch (mode) {
      case 0:
        mat = TransMatrix::Translation(step*i, 0., 0.);
        break;
      case 1:
        mat = TransMatrix::Translation(0., step*i, 0.);
        break;
      case 2:
        mat = TransMatrix::Translation(0., 0., step*i);
        break;
      case 3: {
        //TransMatrix m1 = TransMatrix::Translation(-128,-128,-128);
        //TransMatrix m2 = TransMatrix::Translation(128,128,128);
        TransMatrix m3 = TransMatrix::RotationX(step*i);
        //m3 = m3.MultMatrixRightBy(m1);
        mat = m3; //.MultMatrixLeftBy(m2);
              }
        break;
      case 4:{
        // TransMatrix m1 = TransMatrix::Translation(-128,-128,-128);
        // TransMatrix m2 = TransMatrix::Translation(128,128,128);
        TransMatrix m3 = TransMatrix::RotationY(step*i);
        // m3 = m3.MultMatrixRightBy(m1);
        mat = m3; //.MultMatrixLeftBy(m2);
              }
        break;
      case 5:{
        //TransMatrix m1 = TransMatrix::Translation(-128,-128,-128);
        // TransMatrix m2 = TransMatrix::Translation(128,128,128);
        TransMatrix m3 = TransMatrix::RotationZ(step*i);
        // m3 = m3.MultMatrixRightBy(m1);
        mat = m3; //.MultMatrixLeftBy(m2);
              }
        break;
      case 6: {
                TransMatrix m1 = TransMatrix::Translation(-128,-128,-128);
                TransMatrix m2 = TransMatrix::Translation(128,128,128);
                TransMatrix m3 = TransMatrix::RotationX(step*i);
                TransMatrix m4 = TransMatrix::RotationY(step*i);
                m3 = m3.MultMatrixRightBy(m4);
                m3 = m3.MultMatrixRightBy(m1);
                mat = m3.MultMatrixLeftBy(m2);
              }
              break;
      case 7: {
                TransMatrix m1 = TransMatrix::Translation(-128,-128,-128);
                TransMatrix m2 = TransMatrix::Translation(128,128,128);
                TransMatrix m3 = TransMatrix::RotationX(step*i);
                TransMatrix m4 = TransMatrix::RotationZ(step*i);
                m3 = m3.MultMatrixRightBy(m4);
                m3 = m3.MultMatrixRightBy(m1);
                mat = m3.MultMatrixLeftBy(m2);
              }
              break;
    }

    // Compute the right similarity measure ...
    switch (metric) {
      case 0:
        val = sim->getSSD(mat.GetDataPointer());
        break;
      case 1:
        val = sim->getNCC(mat.GetDataPointer());
        break;
      case 2:
        val = sim->getMI(mat.GetDataPointer());
        break;
      case 3:
        val = sim->getNMI(mat.GetDataPointer());
        break;
    }
    
    std::cout <<  val << ", "; //std::endl;
    // plot it ...
    m_plot->setValue(cnt++, val, useOctCheck->isChecked());
  }
  sim->useOctree(!useOctCheck->isChecked());

  std::cout << " ]" << std::endl << std::endl;

  std::cout << "reg_mi = [ ";
//  m_plot->setPlotRange(iter, step, modesCombo->currentText());
 
  cnt = 0;
  val=0;
  for (int i= -iter; i<iter+1; i++) {
    // Get the right matrix ...
    switch (mode) {
      case 0:
        mat = TransMatrix::Translation(step*i, 0., 0.);
        break;
      case 1:
        mat = TransMatrix::Translation(0., step*i, 0.);
        break;
      case 2:
        mat = TransMatrix::Translation(0., 0., step*i);
        break;
      case 3: {
        // TransMatrix m1 = TransMatrix::Translation(-128,-128,-128);
        // TransMatrix m2 = TransMatrix::Translation(128,128,128);
        TransMatrix m3 = TransMatrix::RotationX(step*i);
        // m3 = m3.MultMatrixRightBy(m1);
        mat = m3; //.MultMatrixLeftBy(m2);
              }
        break;
      case 4:{
        // TransMatrix m1 = TransMatrix::Translation(-128,-128,-128);
        // TransMatrix m2 = TransMatrix::Translation(128,128,128);
        TransMatrix m3 = TransMatrix::RotationY(step*i);
        // m3 = m3.MultMatrixRightBy(m1);
        mat = m3; //.MultMatrixLeftBy(m2);
              }
        break;
      case 5:{
        // TransMatrix m1 = TransMatrix::Translation(-128,-128,-128);
        // TransMatrix m2 = TransMatrix::Translation(128,128,128);
        TransMatrix m3 = TransMatrix::RotationZ(step*i);
        // m3 = m3.MultMatrixRightBy(m1);
        mat = m3; //.MultMatrixLeftBy(m2);
              }
        break;
      case 6: {
                TransMatrix m1 = TransMatrix::Translation(-128,-128,-128);
                TransMatrix m2 = TransMatrix::Translation(128,128,128);
                TransMatrix m3 = TransMatrix::RotationX(step*i);
                TransMatrix m4 = TransMatrix::RotationY(step*i);
                m3 = m3.MultMatrixRightBy(m4);
                m3 = m3.MultMatrixRightBy(m1);
                mat = m3.MultMatrixLeftBy(m2);
              }
              break;
      case 7: {
                TransMatrix m1 = TransMatrix::Translation(-128,-128,-128);
                TransMatrix m2 = TransMatrix::Translation(128,128,128);
                TransMatrix m3 = TransMatrix::RotationX(step*i);
                TransMatrix m4 = TransMatrix::RotationZ(step*i);
                m3 = m3.MultMatrixRightBy(m4);
                m3 = m3.MultMatrixRightBy(m1);
                mat = m3.MultMatrixLeftBy(m2);
              }
              break;

    }

    // Compute the right similarity measure ...
    switch (metric) {
      case 0:
        val = sim->getSSD(mat.GetDataPointer());
        break;
      case 1:
        val = sim->getNCC(mat.GetDataPointer());
        break;
      case 2:
        val = sim->getMI(mat.GetDataPointer());
        break;
      case 3:
        val = sim->getNMI(mat.GetDataPointer());
        break;
    }
    
    std::cout << val << ", "; // std::endl;
    // plot it ...
    m_plot->setValue(cnt++, val, !useOctCheck->isChecked());
  }
  std::cout << " ]" << std::endl << std::endl;

}

void MainWindow::onSave() {
#if QT_VERSION >= 0x040100

#ifndef QT_NO_FILEDIALOG
  const QString fileName = QFileDialog::getSaveFileName(
      this, "Save File Name", QString(),
      "PDF Documents (*.pdf)");
#else
  const QString fileName = "omi.pdf";
#endif

  if ( !fileName.isEmpty() )
  {
    QPrinter printer;
    printer.setOutputFormat(QPrinter::PdfFormat);
    printer.setOrientation(QPrinter::Landscape);
    printer.setOutputFileName(fileName);

    printer.setCreator("Octree Similarity");
    m_plot->print(printer);
  }
#endif
}

