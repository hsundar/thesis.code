#ifndef __PATIENT_DATA_H_
#define __PATIENT_DATA_H_

// Includes ...
#include <QList>
#include <QString>


// forward declarations ...
class Volume;
template<typename T, int n> class Field;


class PatientData {
  public:
    PatientData();
    virtual ~PatientData();

    void add


  protected:  
    QString			m_strName;

    // can have multiple volumes per patient. For now we assume that if there is
    // one volume then it is a single 3D volume, else we have a sequence. We
    // currently do not deal with the case when we have multiple scans of the
    // same patient.
    QList <Volume *>		m_Volume;
    // Can have a set of corresponding deformation fields too ...
    QList <Field<float, 3> *> 	m_defFields;

    // Only one functional Anisotropy/ Principal Direction or diffetion tensor
    // data can be loaded per patient.
    Field <float, 3>		*m_PrincipalDirection;
    Field <float, 1>		*m_FractionalAnisotropy;
    Field <float, 6>		*m_DiffusionTensor;

    // For now the octrees are going to be stroed within the volume.
}

#endif

