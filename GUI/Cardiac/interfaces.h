#ifndef INTERFACES_H
#define INTERFACES_H

class Volume;
class Transformation;

class QString;
class QStringList;

class RegInterface
{
public:
    virtual ~RegInterface() {}

    virtual QStringList algos() const = 0;
    
    virtual void setSubjectImage(Volume *sub) = 0;
    virtual void setTemplateImage(Volume *sub) = 0;
    virtual void setInitialTransformation(Transformation *trans) =0;
    
    virtual void Register(const QString &algo) = 0;

    virtual void getFinalTransformation(Transformation *trans) =0;
};

class FilterInterface
{
public:
    virtual ~FilterInterface() {}

    virtual QStringList filters() const = 0;
    virtual int filterImage(const QString &filter, const Volume *in_vol, Volume *out_vol,
                               QWidget *parent) = 0;
};

Q_DECLARE_INTERFACE(RegInterface,
                    "edu.upenn.sbia.cardiac.RegInterface/1.0")
Q_DECLARE_INTERFACE(FilterInterface,
                    "edu.upenn.sbia.cardiac.FilterInterface/1.0")

#endif
