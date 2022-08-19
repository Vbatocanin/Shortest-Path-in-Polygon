#ifndef GA08_DELAUNAY_TRIANGULATION_H
#define GA08_DELAUNAY_TRIANGULATION_H

#include "algoritambaza.h"
#include "pomocnefunkcije.h"

#include "ga06_dcel.h"

class Delone: public AlgoritamBaza{
public:
    Delone(QWidget *pCrtanje,
                int pauzaKoraka,
                const bool &naivni = false,
                std::string imeDatoteke = "",
                int brojTacaka = BROJ_SLUCAJNIH_OBJEKATA);

    void pokreniAlgoritam() final;
    void crtajAlgoritam(QPainter *painter) const final;
    void pokreniNaivniAlgoritam() final;
    void crtajNaivniAlgoritam(QPainter *painter) const final;

private:
    std::vector<Vertex*> generisiNasumicneTacke(int brojTacaka = BROJ_SLUCAJNIH_OBJEKATA) const;
    std::vector<Vertex*> ucitajPodatkeIzDatoteke(std::string imeDatoteke) const;

    std::vector<Vertex*> _tacke;

};

#endif // GA08_DELAUNAY_TRIANGULATION_H
