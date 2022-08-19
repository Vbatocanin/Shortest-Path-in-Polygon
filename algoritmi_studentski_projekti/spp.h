#ifndef SPP_H
#define SPP_H

#include "algoritambaza.h"
#include "./algoritmi_sa_vezbi/ga07_triangulation.h"
#include <queue>

using PathSet = std::set<Vertex>;

/*Indeksi koji se koriste pri hesiranju pocetne i krajnje tacke u mapi*/
enum ReservedIndex
{
    START_POSITION=-2,
    END_POSITION=-3
};

class SPP : public AlgoritamBaza
{
public:
    /* Konstruktor i destruktor klase */
    SPP(QWidget *,
        int,
        const bool & = false,
        std::string = "",
        int = BROJ_SLUCAJNIH_OBJEKATA);
    virtual ~SPP() override;

    /* Virtuelni metodi iz natklase */
    void pokreniAlgoritam() final;
    void crtajAlgoritam(QPainter *) const final;
    void pokreniNaivniAlgoritam() final;
    void crtajNaivniAlgoritam(QPainter *) const final;

    /* Dodatni metod za grubu silu */
    void pokreniAlgoritamGrubeSile();

    /* Dohvataci rezultata algoritama */
    PathSet getGlavni() const;
    PathSet getNaivni() const;
    PathSet getGruba() const;

    std::vector<QPointF> ucitajNasumicneTacke(int brojTacaka) const;

    std::vector<QPointF> generisiNasumicneTacke(int brojTacaka) const;
    bool checkDiagonal(Vertex *v, Vertex *u);

    const std::unordered_map<Vertex, std::vector<Vertex> > &visibilityGraph() const;
    void setVisibilityGraph(const std::unordered_map<Vertex, std::vector<Vertex> > &newVisibilityGraph);

    std::vector<int> dijkstra(int src);
    int minDistance(std::vector<int> dist, std::vector<bool> sptSet, int V);
    double vDistance(int u, int v);
    bool seceIvicuPoligona(Vertex *u, Vertex *v);
    void triangulatePolygon();
    void generateVisibilityGraphFromVertex(Vertex* v);
    int getAngle(Vertex *a, Vertex *b, Vertex *c);
    void findVertexThatSees(Vertex *start, Vertex *goal);
    bool imaPresek(Vertex *u, Vertex *v);
    bool belongToSameFOV(int ix1, int ix2);
    void tightenImperfectPath();
    Vertex* generateRandomPositionInPolygon();
    bool isVisibleFromTo(Vertex *v, Vertex *u);
    double chainLen(int start, int target);
    std::vector<int> traverseHull(Vertex *start, Vertex *target, bool ccw);
    bool isFasterThan(std::vector<int> path1, std::vector<int> path2);
signals:
    void visibilityGraphChanged();

private:
    std::vector<QPointF> ucitajPodatkeIzDatoteke(std::string imeDatoteke) const;
    DCEL _polygon;
    DCEL _naivePolygon;
    std::vector<std::pair<Vertex *, Vertex *>> _triDiagonals;
    //Pocetna i krajnja tacka
    Vertex* _p;
    Vertex* _q;
    Vertex _pObj;
    Vertex _qObj;
    Vertex* _pNaive;
    Vertex* _qNaive;

    //A visibility graph represented by a map of polygon vertices, where every vertex
    //is associated with a vector of vertices which are visible from their perspective
    std::map<int, std::set<int>> _visibilityGraph;
    void generateVisibilityGraph();
    //Variable that stores the shortest path from p to q
    std::vector<int>  _shortestPathInPolygon;
    std::vector<std::vector<double>> _adjMatrix;

    std::priority_queue<std::pair<int,Vertex*>,
                        std::vector<std::pair<int,Vertex*>>,
                        std::greater<std::pair<int,Vertex*>>> _eventQueue;
    std::priority_queue<std::tuple<double,Vertex*,Vertex*>,
                        std::vector<std::tuple<double,Vertex*,Vertex*>>,
                        std::greater<std::tuple<double,Vertex*,Vertex*>>> _activeEdges;
    bool _monotone;
};

#endif // SPP_H
