#include "qte.h"
#include "implicits.h"
#include "ui_interface.h"

MainWindow::MainWindow() : QMainWindow(), uiw(new Ui::Assets)
{
	// Chargement de l'interface
    uiw->setupUi(this);

	// Chargement du GLWidget
	meshWidget = new MeshWidget;
	QGridLayout* GLlayout = new QGridLayout;
	GLlayout->addWidget(meshWidget, 0, 0);
	GLlayout->setContentsMargins(0, 0, 0, 0);
    uiw->widget_GL->setLayout(GLlayout);

	// Creation des connect
	CreateActions();

	meshWidget->SetCamera(Camera(Vector(10, 0, 0), Vector(0.0, 0.0, 0.0)));
    tx=0;
    ty=0;
    tz=0;

    hx=1.0;
    hy=1.0;
    hz=1.0;

    rx=0;
    ry=0;
    rz=0;

}

MainWindow::~MainWindow()
{
	delete meshWidget;
}

void MainWindow::CreateActions()
{
	// Buttons
    connect(uiw->boxMesh, SIGNAL(clicked()), this, SLOT(BoxMeshExample()));
    connect(uiw->disqueButton, SIGNAL(clicked()), this, SLOT(DisqueMeshExample()));
    connect(uiw->cylindreButton, SIGNAL(clicked()), this, SLOT(CylindreMeshExample()));
    connect(uiw->sphereButton, SIGNAL(clicked()), this, SLOT(SphereMeshExample()));
    connect(uiw->toreButton, SIGNAL(clicked()), this, SLOT(ToreMeshExample()));
    connect(uiw->capsuleButton, SIGNAL(clicked()), this, SLOT(CapsuleMeshExample()));
    connect(uiw->lampeButton, SIGNAL(clicked()), this, SLOT(LampeMeshExample()));
    connect(uiw->sphereImplicit, SIGNAL(clicked()), this, SLOT(SphereImplicitExample()));
    connect(uiw->couronneButton, SIGNAL(clicked()), this, SLOT(MargueriteMeshExample()));
    connect(uiw->resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));
    connect(uiw->wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->deform, SIGNAL(clicked()), this, SLOT(ObjectDeform()));
    connect(uiw->occlusion, SIGNAL(clicked()), this, SLOT(Occlusion()));
    connect(uiw->capture, SIGNAL(clicked()), this, SLOT(Capture()));



	// Widget edition
	connect(meshWidget, SIGNAL(_signalEditSceneLeft(const Ray&)), this, SLOT(editingSceneLeft(const Ray&)));
	connect(meshWidget, SIGNAL(_signalEditSceneRight(const Ray&)), this, SLOT(editingSceneRight(const Ray&)));

    // Deformation
    connect(uiw->XValueD, &QSpinBox::valueChanged, [this](const int& val)->void{uiw->XSliderD->setValue(val);UpdateGeometry();});
    connect(uiw->XSliderD, &QSlider::valueChanged,[this](const int& val)->void{uiw->XValueD->setValue(val);UpdateGeometry();});
    connect(uiw->YValueD, &QSpinBox::valueChanged, [this](const int& val)->void{uiw->YSliderD->setValue(val);UpdateGeometry();});
    connect(uiw->YSliderD, &QSlider::valueChanged,[this](const int& val)->void{uiw->YValueD->setValue(val);UpdateGeometry();});
    connect(uiw->ZValueD, &QSpinBox::valueChanged, [this](const int& val)->void{uiw->ZSliderD->setValue(val);UpdateGeometry();});
    connect(uiw->ZSliderD, &QSlider::valueChanged,[this](const int& val)->void{uiw->ZValueD->setValue(val);UpdateGeometry();});
    connect(uiw->RValueD, &QDoubleSpinBox::valueChanged, [this](const double& val)->void{uiw->RSliderD->setValue(val*10);UpdateGeometry();});
    connect(uiw->RSliderD, &QSlider::valueChanged,[this](const int& val)->void{uiw->RValueD->setValue(val/10.0);UpdateGeometry();});

    // Translation
    connect(uiw->XValueT, &QDoubleSpinBox::valueChanged, [this](const double& val)->void{uiw->XSliderT->setValue(val*10);TranslationSlider();});
    connect(uiw->XSliderT, &QSlider::valueChanged,[this](const int& val)->void{uiw->XValueT->setValue(val/10.0);TranslationSlider();});
    connect(uiw->YValueT, &QDoubleSpinBox::valueChanged, [this](const double& val)->void{uiw->YSliderT->setValue(val*10);TranslationSlider();});
    connect(uiw->YSliderT, &QSlider::valueChanged,[this](const int& val)->void{uiw->YValueT->setValue(val/10.0);TranslationSlider();});
    connect(uiw->ZValueT, &QDoubleSpinBox::valueChanged, [this](const double& val)->void{uiw->ZSliderT->setValue(val*10);TranslationSlider();});
    connect(uiw->ZSliderT, &QSlider::valueChanged,[this](const int& val)->void{uiw->ZValueT->setValue(val/10.0);TranslationSlider();});
    connect(uiw->resetT, SIGNAL(clicked()), this, SLOT(ResetTranslationSlider()));


    // Homothétie
    connect(uiw->XValueH, &QDoubleSpinBox::valueChanged, [this](const double& val)->void{uiw->XSliderH->setValue(val*100);HomothetieSlider();});
    connect(uiw->XSliderH, &QSlider::valueChanged,[this](const int& val)->void{uiw->XValueH->setValue(val/100.0);HomothetieSlider();});
    connect(uiw->YValueH, &QDoubleSpinBox::valueChanged, [this](const double& val)->void{uiw->YSliderH->setValue(val*100);HomothetieSlider();});
    connect(uiw->YSliderH, &QSlider::valueChanged,[this](const int& val)->void{uiw->YValueH->setValue(val/100.0);HomothetieSlider();});
    connect(uiw->ZValueH, &QDoubleSpinBox::valueChanged, [this](const double& val)->void{uiw->ZSliderH->setValue(val*100);HomothetieSlider();});
    connect(uiw->ZSliderH, &QSlider::valueChanged,[this](const int& val)->void{uiw->ZValueH->setValue(val/100.0);HomothetieSlider();});
    connect(uiw->resetH, SIGNAL(clicked()), this, SLOT(ResetHomothetieSlider()));

    // Rotation
    connect(uiw->XValueR, &QDoubleSpinBox::valueChanged, [this](const double& val)->void{ uiw->XSliderR->setValue(val);RotationSlider();});
    connect(uiw->XSliderR, &QSlider::valueChanged,[this](const int& val)->void{uiw->XValueR->setValue(val);RotationSlider();});
    connect(uiw->YValueR, &QDoubleSpinBox::valueChanged, [this](const double& val)->void{uiw->YSliderR->setValue(val);RotationSlider();});
    connect(uiw->YSliderR, &QSlider::valueChanged,[this](const int& val)->void{uiw->YValueR->setValue(val);RotationSlider();});
    connect(uiw->ZValueR, &QDoubleSpinBox::valueChanged, [this](const double& val)->void{uiw->ZSliderR->setValue(val);RotationSlider();});
    connect(uiw->ZSliderR, &QSlider::valueChanged,[this](const int& val)->void{uiw->ZValueR->setValue(val);RotationSlider();});
    connect(uiw->resetR, SIGNAL(clicked()), this, SLOT(ResetRotationSlider()));

}

void MainWindow::editingSceneLeft(const Ray&)
{
}

void MainWindow::editingSceneRight(const Ray&)
{
}

void MainWindow::BoxMeshExample()
{
    ResetAll();
	Mesh boxMesh = Mesh(Box(1.0));

	std::vector<Color> cols;
	cols.resize(boxMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

	meshColor = MeshColor(boxMesh, cols, boxMesh.VertexIndexes());
    ResetAll();
	UpdateGeometry();
}

/*!
 * \brief MainWindow::DisqueMeshExample Création d'un Mesh de disque
 */
void MainWindow::DisqueMeshExample()
{
    ResetAll();
    int nb= uiw->nbPointDisque->value();
    Mesh disqueMesh = Mesh(Disque(Vector(0,0,0),3),nb);
    meshColor = MeshColor(disqueMesh);
    ResetAll();
    UpdateGeometry();

}

/*!
 * \brief MainWindow::CylindreMeshExample Création d'un Mesh d'un cylindre
 */
void MainWindow::CylindreMeshExample()
{
    ResetAll();
    int nbP= uiw->nbPointCylindre->value();
    int nbC= uiw->nbCercleCylindre->value();
    Mesh cylindreMesh = Mesh(Cylindre(Vector(0,0,0),3,4),nbP,nbC,uiw->fermeCyl->isChecked());
    meshColor = MeshColor(cylindreMesh);
    ResetAll();
    UpdateGeometry();

}

/*!
 * \brief MainWindow::SphereMeshExample Création d'un Mesh d'une sphère
 */
void MainWindow::SphereMeshExample()
{
    ResetAll();
    int nbP= uiw->nbPointSphere->value();
    int nbC= uiw->nbCercleSphere->value();
    Mesh sphereMesh = Mesh(Sphere(),nbP,nbC,uiw->demiSphere->isChecked());
    meshColor = MeshColor(sphereMesh);
    UpdateGeometry();

}

/*!
 * \brief MainWindow::ToreMeshExample Création d'un Mesh d'un tore
 */
void MainWindow::ToreMeshExample()
{
    ResetAll();
    int nbP= uiw->nbPointTore->value();
    int nbC= uiw->nbCercleTore->value();
    Mesh toreMesh = Mesh(Tore(Vector(0,0,0),3,1),nbP,nbC);
    meshColor = MeshColor(toreMesh);
    ResetAll();
    UpdateGeometry();

}

/*!
 * \brief MainWindow::CapsuleMeshExample Création d'un Mesh d'une capsule
 */
void MainWindow::CapsuleMeshExample()
{
    ResetAll();
    int nbP= uiw->nbPointCaps->value();
    int nbCD= uiw->nbCercleCapsDem->value();
    int nbCC= uiw->nbCercleCapsCyl->value();

    Mesh capsuleMesh = Mesh(Capsule(Vector(0,0,0),2,8),nbP,nbCC,nbCD);

    meshColor = MeshColor(capsuleMesh);
    UpdateGeometry();

}

/*!
 * \brief MainWindow::SphereImplicitExample Création d'un Mesh d'une sphère
 */
void MainWindow::SphereImplicitExample()
{
  ResetAll();
  AnalyticScalarField implicit;

  Mesh implicitMesh;
  implicit.Polygonize(31, implicitMesh, Box(2.0));

  std::vector<Color> cols;
  cols.resize(implicitMesh.Vertexes());
  for (size_t i = 0; i < cols.size(); i++)
    cols[i] = Color(0.8, 0.8, 0.8);

  meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
  ResetAll();
  UpdateGeometry();
}

/*!
 * \brief MainWindow::LampeMeshExample Création d'un Mesh d'une lampe
 */
void MainWindow::LampeMeshExample()
{
    ResetAll();
    int nbPS= uiw->nbPointLampSpher->value();
    int nbCS= uiw->nbCercleLampSpher->value();
    int nbPC= uiw->nbPointLampCyl->value();
    int nbCC= uiw->nbCercleLampCyl->value();
    int i;

    Mesh domeMesh = Mesh(Sphere(Vector(0,0,0),5),nbPS,nbCS);
    std::vector<Color> cols;
    for (i = 0; i < domeMesh.Vertexes(); i++)
        cols.push_back(Color(1.0,1.0, 1.0));
    domeMesh.SphereWarp(Sphere(Vector(0,0,-8),7),Vector(0,0,12));
    Mesh piedMesh=Mesh(Cylindre(Vector(0,0,-2.66),0.2,14.75),nbPC,nbCC);
    domeMesh.Merge(piedMesh);
    for (; i < domeMesh.Vertexes(); i++)
        cols.push_back(Color(255,215, 0));
    Mesh socleMesh=Mesh(Cylindre(Vector(0,0,-10.54),2.5,1),nbPC,nbCC);
    domeMesh.Merge(socleMesh);
    for (; i < domeMesh.Vertexes(); i++)
        cols.push_back(Color(255,215, 0));
    meshColor = MeshColor(domeMesh, cols, domeMesh.VertexIndexes());
    ResetAll();
    UpdateGeometry();
}

/*!
 * \brief MainWindow::MargueriteMeshExample Création d'un mesh d'une marguerite superposée d'une boule
 */
void MainWindow::MargueriteMeshExample(){
    ResetAll();

    int nbS = uiw->nbSphereCouronne->value();
    int nbP = uiw->nbPointCouronne->value();
    int nbC = uiw->nbCercleCouronne->value();

    Mesh sphereMesh = Mesh(Sphere(),nbP,nbC);
    Mesh spheretmp;

    Mesh sphereMesh2 = Mesh(Sphere(Vector(0,0,6.5),3),nbP,nbC);
    sphereMesh.Merge(sphereMesh2);

    if(nbS > 0 && nbP > 0 && nbC > 0){

        float step2 = 2.0 * M_PI / (nbS);
        float alpha2;
        Vector c;
        for(int i=0;i<nbS;i++){
            alpha2 = i*step2;
            c = Vector(cos(alpha2),sin(alpha2),0)*3;
            spheretmp = Mesh(Sphere(c,2),nbP,nbC);
            sphereMesh.Merge(spheretmp);

        }

        meshColor = MeshColor(sphereMesh);
        UpdateGeometry();
    }

}

void MainWindow::UpdateGeometry()
{
	meshWidget->ClearAll();
    meshWidget->AddMesh("BoxMesh", meshColor);
    if(uiw->RValueD->value() > 0){
        meshWidget->AddMesh("Sphere",Mesh(Sphere(Vector(uiw->XValueD->value(),uiw->YValueD->value(),uiw->ZValueD->value()),uiw->RValueD->value()),10,10));
    }

    uiw->lineEdit->setText(QString::number(meshColor.Vertexes()));
    uiw->lineEdit_2->setText(QString::number(meshColor.Triangles()));

	UpdateMaterial();
}

void MainWindow::TranslationSlider(){
    meshColor.transfer(Vector(uiw->XValueT->value()-tx,uiw->YValueT->value()-ty,uiw->ZValueT->value()-tz));
    tx = uiw->XValueT->value();
    ty = uiw->YValueT->value();
    tz = uiw->ZValueT->value();
    UpdateGeometry();
}

void MainWindow::HomothetieSlider(){
    meshColor.transform(Matrix::MatHomothetie(uiw->XValueH->value()/hx,uiw->YValueH->value()/hy,uiw->ZValueH->value()/hz));

    hx = uiw->XValueH->value();
    hy = uiw->YValueH->value();
    hz = uiw->ZValueH->value();
    UpdateGeometry();
}

void MainWindow::RotationSlider(){
    meshColor.transform(Matrix::MatRotX(Math::DegreeToRadian(uiw->XValueR->value()-rx)));
    meshColor.transform(Matrix::MatRotY(Math::DegreeToRadian(uiw->YValueR->value()-ry)));
    meshColor.transform(Matrix::MatRotZ(Math::DegreeToRadian(uiw->ZValueR->value()-rz)));

    rx = uiw->XValueR->value();
    ry = uiw->YValueR->value();
    rz = uiw->ZValueR->value();
    UpdateGeometry();
}

void MainWindow::ResetTranslationSlider(){
    uiw->XSliderT->setValue(0);
    uiw->XValueT->setValue(0);

    uiw->YSliderT->setValue(0);
    uiw->YValueT->setValue(0);

    uiw->ZSliderT->setValue(0);
    uiw->ZValueT->setValue(0);

    TranslationSlider();

    tx = 0;
    ty = 0;
    tz = 0;

}

void MainWindow::ResetHomothetieSlider(){
    uiw->XSliderH->setValue(100);
    uiw->XValueH->setValue(1);

    uiw->YSliderH->setValue(100);
    uiw->YValueH->setValue(1);

    uiw->ZSliderH->setValue(100);
    uiw->ZValueH->setValue(1);

    HomothetieSlider();

    hx = 1;
    hy = 1;
    hz = 1;

}

void MainWindow::ResetRotationSlider(){
    uiw->XSliderR->setValue(0);
    uiw->XValueR->setValue(0);

    uiw->YSliderR->setValue(0);
    uiw->YValueR->setValue(0);

    uiw->ZSliderR->setValue(0);
    uiw->ZValueR->setValue(0);

    RotationSlider();

    rx = 0;
    ry = 0;
    rz = 0;

}

void MainWindow::ResetAll(){
    ResetTranslationSlider();
    ResetHomothetieSlider();
    ResetRotationSlider();
}


void MainWindow::UpdateMaterial()
{
    meshWidget->UseWireframeGlobal(uiw->wireframe->isChecked());

    if (uiw->radioShadingButton_1->isChecked())
		meshWidget->SetMaterialGlobal(MeshMaterial::Normal);
	else
        meshWidget->SetMaterialGlobal(MeshMaterial::Color);

    meshWidget->SetMaterial("Sphere",MeshMaterial::Color);

}

void MainWindow::ResetCamera()
{
	meshWidget->SetCamera(Camera(Vector(-10.0), Vector(0.0)));
}

void MainWindow::ObjectDeform(){
    meshColor.SphereWarp(Sphere(Vector(uiw->XValueD->value(),uiw->YValueD->value(),uiw->ZValueD->value()),uiw->RValueD->value()),Vector(uiw->XVector->value(),uiw->YVector->value(),uiw->ZVector->value()));
    UpdateGeometry();
}

void MainWindow::Occlusion(){
    meshColor.CalculOcclusion(uiw->nbRayon->value());
    uiw->radioShadingButton_1->setChecked(false);
    uiw->radioShadingButton_2->setChecked(true);
    UpdateGeometry();
}

void MainWindow::Capture(){
    meshWidget->SaveScreen();
}
