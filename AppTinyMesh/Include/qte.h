#ifndef __Qte__
#define __Qte__

#include <QtWidgets/qmainwindow.h>
#include "realtime.h"
#include "meshcolor.h"

QT_BEGIN_NAMESPACE
	namespace Ui { class Assets; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
  Q_OBJECT
private:
  Ui::Assets* uiw;           //!< Interface

  MeshWidget* meshWidget;   //!< Viewer
  MeshColor meshColor;		//!< Mesh.

  double tx,ty,tz;
  double hx,hy,hz;
  double rx,ry,rz;

public:
  MainWindow();
  ~MainWindow();
  void CreateActions();
  void UpdateGeometry();
  void DeformationSlider();
  void TranslationSlider();
  void HomothetieSlider();
  void RotationSlider();


public slots:
  void editingSceneLeft(const Ray&);
  void editingSceneRight(const Ray&);
  void BoxMeshExample();
  void DisqueMeshExample();
  void CylindreMeshExample();
  void SphereMeshExample();
  void ToreMeshExample();
  void CapsuleMeshExample();
  void SphereImplicitExample();
  void LampeMeshExample();
  void MargueriteMeshExample();
  void ResetCamera();
  void UpdateMaterial();
  void ObjectDeform();
  void Occlusion();

  void ResetTranslationSlider();
  void ResetHomothetieSlider();
  void ResetRotationSlider();
  void ResetAll();
  void Capture();
};

#endif
