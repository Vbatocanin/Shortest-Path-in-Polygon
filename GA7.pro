QT       += core gui opengl charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
LIBS += -lglut -lGLU
CONFIG += c++14 console

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    algoritambaza.cpp \
    algoritmi_sa_vezbi/ga00_demoiscrtavanja.cpp \
    algoritmi_sa_vezbi/ga01_brisucaprava.cpp \
    algoritmi_sa_vezbi/ga02_3discrtavanje.cpp \
    algoritmi_sa_vezbi/ga03_konveksniomotac.cpp \
    algoritmi_sa_vezbi/ga04_konveksniomotac3d.cpp \
    algoritmi_sa_vezbi/ga05_preseciduzi.cpp \
    algoritmi_sa_vezbi/ga06_dcel.cpp \
    algoritmi_sa_vezbi/ga06_dceldemo.cpp \
    algoritmi_sa_vezbi/ga07_triangulation.cpp \
    algoritmi_sa_vezbi/ga08_delaunay_trinagulation.cpp \
    algoritmi_studentski_projekti/ga06_presekPravougaonika.cpp \
    algoritmi_studentski_projekti/spp.cpp \
    animacijanit.cpp \
    main.cpp \
    mainwindow.cpp \
    oblastcrtanja.cpp \
    oblastcrtanjaopengl.cpp \
    pomocnefunkcije.cpp \
    timemeasurementthread.cpp

HEADERS += \
    algoritambaza.h \
    algoritmi_sa_vezbi/ga00_demoiscrtavanja.h \
    algoritmi_sa_vezbi/ga01_brisucaprava.h \
    algoritmi_sa_vezbi/ga02_3discrtavanje.h \
    algoritmi_sa_vezbi/ga03_konveksniomotac.h \
    algoritmi_sa_vezbi/ga04_konveksni3dDatastructures.h \
    algoritmi_sa_vezbi/ga04_konveksniomotac3d.h \
    algoritmi_sa_vezbi/ga05_datastructures.h \
    algoritmi_sa_vezbi/ga05_preseciduzi.h \
    algoritmi_sa_vezbi/ga06_dcel.h \
    algoritmi_sa_vezbi/ga06_dceldemo.h \
    algoritmi_sa_vezbi/ga07_datastructures.h \
    algoritmi_sa_vezbi/ga07_triangulation.h \
    algoritmi_sa_vezbi/ga08_delaunay_triangulation.h \
    algoritmi_studentski_projekti/ga06_presekPravougaonika.h \
    algoritmi_studentski_projekti/spp.h \
    animacijanit.h \
    config.h \
    mainwindow.h \
    oblastcrtanja.h \
    oblastcrtanjaopengl.h \
    pomocnefunkcije.h \
    timemeasurementthread.h

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

DISTFILES += \
    docs/ga06_presekPravougaonika.docx \
    docs/ga06_presekPravougaonika.pdf \
    input_files/ga00_DCELDemo/mushroom.off \
    input_files/ga00_DCELDemo/test0.off \
    input_files/ga00_DCELDemo/test1.off \
    input_files/ga00_KonveksniOmotac3D/input1.txt \
    input_files/ga00_PreseciDuzi/input1.txt \
    input_files/ga00_Triangulation/input1.txt \
    input_files/ga00_Triangulation/input2.txt \
    input_files/ga06_PresekPravougaonika/input1.txt \
    input_files/ga06_PresekPravougaonika/input2.txt \
    input_files/ga06_PresekPravougaonika/input3.txt
