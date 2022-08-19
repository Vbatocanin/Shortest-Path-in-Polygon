#include "spp.h"
#include <fstream>

SPP::SPP(QWidget *pCrtanje,
         int pauzaKoraka,
         const bool &naivni,
         std::string imeDatoteke,
         int brojTacaka)
    :AlgoritamBaza(pCrtanje, pauzaKoraka, naivni)
{
    std::vector<QPointF> tacke;

    /* Ucitavanje tacka */
    if (imeDatoteke != "")
        tacke = ucitajPodatkeIzDatoteke(imeDatoteke);
    else
        tacke = ucitajNasumicneTacke(brojTacaka);


    _polygon.loadData(tacke);
    _naivePolygon.loadData(tacke);

    for(auto i = 0; i < _polygon.vsize(); i++){
        _polygon.vertex(i)->setId(i);
    }

    _p = generateRandomPositionInPolygon();
    _q = generateRandomPositionInPolygon();

    _pNaive = _p;
    _qNaive = _q;

    for(auto i = 0; i < _naivePolygon.vsize(); i++){
        _naivePolygon.vertex(i)->setId(i);
    }


}

SPP::~SPP(){}

void SPP::pokreniAlgoritam()
{
    int V = _polygon.vsize();

    /* O(V): Setting the Id for all vertices in _polygon*/
    for(auto i = 0; i < V; i++){
        _polygon.vertex(i)->setId(i);

        //Adding vertex i's two adjacent vertices to the adj matrix to later use for dijkstra
        int v1 = i+1 == V ? 0 : i+1;
        int v2 = i-1 == -1 ? V-1 : i-1;

        _visibilityGraph[i].insert(v1);
        _visibilityGraph[i].insert(v2);

        _adjMatrix.emplace_back(std::vector<double>(V));
        for (int j=0; j<V; j++)
        {
            if(i==j)
                _adjMatrix[i][j] = 0;
            else
                _adjMatrix[i][j] = MAXFLOAT;
        }
    }

    if (!isVisibleFromTo(_p,_q))
    {
        auto it = _visibilityGraph[_p->getId()].cbegin();
        auto pathccw = traverseHull(_p, _q, true);
        auto pathcw = traverseHull(_p, _q, false);

        // Cleaning up clockwise and counterclockwise traversals of our convex hull
        pathccw.emplace(pathccw.begin(),_p->getId());
        pathcw.emplace(pathcw.begin(),_p->getId());

        auto itStart = std::find(pathccw.begin()+1, pathccw.end(), _p->getId());
        if(itStart != pathccw.end())
        {
            pathccw = {itStart, pathccw.end()};
        }

        itStart = std::find(pathcw.begin()+1, pathcw.end(), _p->getId());
        if(itStart != pathcw.end())
        {
            pathcw = {itStart, pathcw.end()};
        }

        // Trying out the ccw traversal first
        for (auto vtx: pathccw)
        {
            _shortestPathInPolygon.emplace_back(vtx);
            AlgoritamBaza_updateCanvasAndBlock()
        }
        tightenImperfectPath();
        AlgoritamBaza_updateCanvasAndBlock();

        // Trying out the cw traversal next
        pathccw = _shortestPathInPolygon;
        _shortestPathInPolygon.clear();
        for (auto vtx: pathcw)
        {
            _shortestPathInPolygon.emplace_back(vtx);
            AlgoritamBaza_updateCanvasAndBlock()
        }
        tightenImperfectPath();
        AlgoritamBaza_updateCanvasAndBlock();

        pathcw = _shortestPathInPolygon;
        _shortestPathInPolygon.clear();

        // Choosing better traversal
        if(isFasterThan(pathccw, pathcw))
        {
            for (auto vtx: pathccw)
            {
                _shortestPathInPolygon.emplace_back(vtx);
                AlgoritamBaza_updateCanvasAndBlock();
            }
        }
        else
        {
            for (auto vtx: pathcw)
            {
                _shortestPathInPolygon.emplace_back(vtx);
                AlgoritamBaza_updateCanvasAndBlock();
            }
        }
        AlgoritamBaza_updateCanvasAndBlock();
    }
    else
    {
        _shortestPathInPolygon = {_p->getId(),_q->getId()};
        AlgoritamBaza_updateCanvasAndBlock()
    }

}

bool SPP::isFasterThan(std::vector<int> path1, std::vector<int> path2)
{
    double distance1 = 0;
    int path1Len = path1.size();

    double distance2 = 0;
    int path2Len = path2.size();

    for(int i=0; i<path1Len-1; ++i)
    {
        distance1 += vDistance(path1[i],path1[i+1]);
    }
    for(int i=0; i<path2Len-1; ++i)
    {
        distance2 += vDistance(path2[i],path2[i+1]);
    }

    return distance1<distance2;
}

void SPP::tightenImperfectPath()
{
    std::vector<int> newPath;

    Vertex* currentVtx1;
    Vertex* currentVtx2;
    int oldLen = _shortestPathInPolygon.size();
    std::vector<int> oldPath = _shortestPathInPolygon;
    _shortestPathInPolygon.clear();


    if(oldLen <= 2)
    {
        _shortestPathInPolygon = oldPath;
        return;
    }
    else if(oldLen == 3)
    {
        auto vtx1 = _polygon.vertex(oldPath[0]);
        auto vtx3 = _polygon.vertex(oldPath[2]);
        if(!isVisibleFromTo(vtx1,vtx3))
        {
            _shortestPathInPolygon = oldPath;
        }
        else
        {
            _shortestPathInPolygon.emplace_back(vtx1->getId());
            _shortestPathInPolygon.emplace_back(vtx3->getId());
        }
        return;
    }

    int i=0;
    int j=oldLen-1;
    bool pathAltered = false;
    for(;;)
    {
        if (i>=oldLen-2)break;
        currentVtx1 = _polygon.vertex(oldPath[i]);
        pathAltered = false;

        for(;;)
        {

            if(i>=j-1)
            {
                if(_shortestPathInPolygon.back() != currentVtx1->getId())
                    _shortestPathInPolygon.emplace_back(currentVtx1->getId());
                break;
            }
            currentVtx2 = _polygon.vertex(oldPath[j]);

            if(isVisibleFromTo(currentVtx1,currentVtx2))
            {
                if(_shortestPathInPolygon.back() != currentVtx1->getId())
                    _shortestPathInPolygon.emplace_back(currentVtx1->getId());
                _shortestPathInPolygon.emplace_back(currentVtx2->getId());
                i = j;
                j=oldLen-1;

                currentVtx1 = currentVtx2;
                pathAltered = true;
                AlgoritamBaza_updateCanvasAndBlock();

                continue;
            }
            --j;
        }
        AlgoritamBaza_updateCanvasAndBlock();

        ++i;
        j=oldLen-1;

        if(pathAltered == false && _shortestPathInPolygon.back() != currentVtx1->getId())
            _shortestPathInPolygon.emplace_back(currentVtx1->getId());
    }



    if(_shortestPathInPolygon.back() != _q->getId())
        _shortestPathInPolygon.emplace_back(_q->getId());
    AlgoritamBaza_updateCanvasAndBlock();
}

std::vector<int> SPP::traverseHull(Vertex* start, Vertex* target, bool ccw)
{

    std::vector<int> path1;
    std::vector<int> path2;

    Vertex* currentVtx = start;

    // Counterclockwise
    if(ccw)
    {
        for (;;)
        {
            path1.emplace_back(currentVtx->getId());

            if (currentVtx->getId() == target->getId())break;
            currentVtx = currentVtx->incidentEdge()->next()->origin();
        }

        return path1;
    }

    // Clockwise
    for (;;)
    {
        path2.emplace_back(currentVtx->getId());

        if (currentVtx == target)break;
        currentVtx = currentVtx->incidentEdge()->prev()->origin();
    }

    return path2;
}

bool SPP::isVisibleFromTo(Vertex* v, Vertex* u)
{
    if(_visibilityGraph[v->getId()].find(u->getId()) != _visibilityGraph[v->getId()].end())
        return true;

    QPointF presek;

    /* Provera da li je (v, u) spoljasnja dijagonala */
    bool badDiag = checkDiagonal(v, u);

    /* Ako jeste, duz (v, u) nije odgovarajuca, prelazi se na sledecu duz */
    if (badDiag) return false;

    /* Provera da li dijagonala sece neku od ivica poligona */
    for(auto k = 0ul; k < _polygon.esize()/2; k++){
        auto edge = _polygon.edge(k);
        if (edge == v->incidentEdge() || edge == v->incidentEdge()->prev() ||
            edge == u->incidentEdge() || edge == u->incidentEdge()->prev())
            continue;

        if (pomocneFunkcije::presekDuzi(QLineF(edge->origin()->coordinates(),
                                               edge->twin()->origin()->coordinates()),
                                        QLineF(v->coordinates(), u->coordinates()),
                                        presek)){
            /* Postoji presek */
            badDiag = true;
            break;
        }
    }

    /* Ako sece, duz (v, u) nije odgovarajuca, prelazi se na sledecu duz */
    if (badDiag) return false;

    /* Ako ispunjava sve uslove, dodajemo dijagonale (v, u) i (u, v) u graf vidljivosti*/
    if (!badDiag){
        _visibilityGraph[v->getId()].insert(u->getId());
        _visibilityGraph[u->getId()].insert(v->getId());
    }


    if (_adjMatrix[v->getId()][u->getId()]==MAXFLOAT)
    {
        double dist = vDistance(v->getId(),u->getId());
        _adjMatrix[v->getId()][u->getId()] = dist;
        _adjMatrix[u->getId()][v->getId()] = dist;
    }
    return true;
}


void SPP::generateVisibilityGraphFromVertex(Vertex* v)
{
    /* Slozenost algoritma: O(n^2) */
    auto i = v->getId();
    for (auto j = 0; j < _polygon.vsize(); j++){
        if (i == j)continue;

        auto u = _polygon.vertex(j);
        if(_visibilityGraph[v->getId()].find(u->getId())!=_visibilityGraph[v->getId()].end())continue;

        QPointF presek;

        /* Provera da li je (v, u) spoljasnja dijagonala */
        bool badDiag = checkDiagonal(v, u);

        /* Ako jeste, duz (v, u) nije odgovarajuca, prelazi se na sledecu duz */
        if (badDiag) continue;

        /* Provera da li dijagonala sece neku od ivica poligona */
        for(auto k = 0ul; k < _polygon.esize()/2; k++){
            auto edge = _polygon.edge(k);
            if (edge == v->incidentEdge() || edge == v->incidentEdge()->prev() ||
                edge == u->incidentEdge() || edge == u->incidentEdge()->prev())
                continue;

            if (pomocneFunkcije::presekDuzi(QLineF(edge->origin()->coordinates(),
                                                   edge->twin()->origin()->coordinates()),
                                            QLineF(v->coordinates(), u->coordinates()),
                                            presek)){
                /* Postoji presek */
                badDiag = true;
                break;
            }
        }

        /* Ako sece, duz (v, u) nije odgovarajuca, prelazi se na sledecu duz */
        if (badDiag) continue;

        /* Ako ispunjava sve uslove, dodajemo dijagonale (v, u) i (u, v) u graf vidljivosti*/
        if (!badDiag){
            _visibilityGraph[v->getId()].insert(u->getId());
            _visibilityGraph[u->getId()].insert(v->getId());
            AlgoritamBaza_updateCanvasAndBlock()
        }
    }

    _visibilityGraph[v->getId()].insert(v->getId());

    for (auto u: _visibilityGraph[v->getId()])
    {
        if (_adjMatrix[v->getId()][u]!=MAXFLOAT)continue;

        double dist = vDistance(v->getId(),u);
        _adjMatrix[v->getId()][u] = dist;
        _adjMatrix[u][v->getId()] = dist;
    }
    AlgoritamBaza_updateCanvasAndBlock()
}

void SPP::crtajAlgoritam(QPainter *painter) const
{
    if (!painter) return;

    QPen regular = painter->pen();
    regular.setWidth(3);

    QPen green = painter->pen();
    green.setColor(Qt::darkGreen);
    green.setWidth(3);
    painter->setPen(green);

    painter->drawEllipse(_p->coordinates(),20,20);
    painter->drawEllipse(_q->coordinates(),20,20);


    // Drawing the first part, middle part, and finally the end segments for the shortest path from p to q
    int pathSize = _shortestPathInPolygon.size();
    pathSize = pathSize > 0 ? pathSize : 0;

    if(pathSize==0)
    {
        for (auto elem: _visibilityGraph)
        {
            Vertex* vtx = _polygon.vertex(elem.first);

            painter->drawText(vtx->coordinates(),QString::number(vtx->getId()));

            for(auto vertex: elem.second)
            {
                painter->drawLine(_polygon.vertex(elem.first)->coordinates(),
                                  _polygon.vertex(vertex)->coordinates());
            }
        }
    }
    else
    {
        for(int i=0; i<pathSize-1; i++)
        {
            painter->drawText(_polygon.vertex(_shortestPathInPolygon[i])->coordinates(),
                              QString::number(_polygon.vertex(_shortestPathInPolygon[i])->getId()));
            painter->drawLine(_polygon.vertex(_shortestPathInPolygon[i])->coordinates(),
                          _polygon.vertex(_shortestPathInPolygon[i+1])->coordinates());
        }
    }


    painter->setPen(regular);
    for (auto v: _polygon.vertices())
    {
        if(_shortestPathInPolygon.size()>1)
        {
            int id1 = v->getId();
            int id2 = v->incidentEdge()->twin()->origin()->getId();
            std::vector<int> ids1{id1,id2};
            std::vector<int> ids2{id2,id1};

            auto ind1 = std::search(_shortestPathInPolygon.begin(), _shortestPathInPolygon.end(), ids1.begin(), ids1.end()) != _shortestPathInPolygon.end();
            auto ind2 = std::search(_shortestPathInPolygon.begin(), _shortestPathInPolygon.end(), ids2.begin(), ids2.end()) != _shortestPathInPolygon.end();

            if(ind1 || ind2)continue;
        }
        painter->drawLine(v->coordinates(), v->incidentEdge()->twin()->origin()->coordinates());
    }
}

void SPP::pokreniNaivniAlgoritam()
{   
    int V = _naivePolygon.vsize();
    for(int i=0; i<V; i++)
    {
        _adjMatrix.emplace_back(std::vector<double>(V));
        for (int j=0; j<V; j++)
        {
            if(i==j)
                _adjMatrix[i][j] = 0;
            else
                _adjMatrix[i][j] = MAXFLOAT;
        }
    }

    generateVisibilityGraph();

    std::vector<int> currentPath{};

    _shortestPathInPolygon = dijkstra(_pNaive->getId());

    AlgoritamBaza_updateCanvasAndBlock()
    emit animacijaZavrsila();
}

std::vector<int> SPP::dijkstra(int src)
{
    const int V = _naivePolygon.vsize();
    // dist[i] is the shortest
    // distance from src to i
    std::vector<int> dist(V);

    // sptSet[i] will true if vertex i is included / in
    // shortest path tree or shortest distance from src to i
     std::vector<bool> sptSet(V);

    // Parent array to store shortest path tree
    std::vector<int> parent(V);

    for(int i =0; i<V; i++)
    {
        sptSet[i]=false;
        parent[i]=-1;
    }

    // Initialize all distances as INT_MAX
    for (int i = 0; i < V; i++)
        dist[i] = INT_MAX;

    // Distance of source vertex from itself is always 0
    dist[src] = 0;

    // Calculating minimum distances from src to all vertices
    for (int count = 0; count < V - 1; count++) {
        int u = minDistance(dist, sptSet, V);
        sptSet[u] = true;
        for (int v = 0; v < V; v++)
            if (!sptSet[v] && _adjMatrix[u][v]
                && dist[u] + _adjMatrix[u][v] < dist[v]) {
                parent[v] = u;
                dist[v] = dist[u] + _adjMatrix[u][v];
            }
    }

    //Shortest path from q to p reconstruction
    std::vector<int> shortestPath;
    int currentVtx = _qNaive->getId();
    shortestPath.emplace_back(currentVtx);
    while(true)
    {
        currentVtx = parent[currentVtx];
        if(currentVtx == -1)
            break;
        shortestPath.emplace(shortestPath.begin(),currentVtx);
        if(_visibilityGraph[currentVtx].find(_pNaive->getId()) != _visibilityGraph[currentVtx].end())
        {
            shortestPath.emplace(shortestPath.begin(),_pNaive->getId());
            break;
        }
    }
    return shortestPath;
}

// Dijkstra utility function
int SPP::minDistance(std::vector<int> dist, std::vector<bool> sptSet, int V)
{
    // Initialize min value
    int min = INT_MAX, min_index;
    for (int i = 0; i <V; i++)
        if (sptSet[i] == false && dist[i] <= min)
            min = dist[i], min_index = i;
    return min_index;
}


void SPP::crtajNaivniAlgoritam(QPainter *painter) const
{
    if (!painter) return;

    QPen regular = painter->pen();
    regular.setWidth(3);

    QPen green = painter->pen();
    green.setColor(Qt::darkGreen);
    green.setWidth(2);

    painter->setPen(green);

    painter->drawEllipse(_pNaive->coordinates(),20,20);
    painter->drawEllipse(_qNaive->coordinates(),20,20);

    if(_shortestPathInPolygon.size()>1)
    {
        painter->setPen(green);
        // Drawing the first part, middle part, and finally the end segments for the shortest path from p to q
        painter->drawLine(_pNaive->coordinates(), _naivePolygon.vertex(_shortestPathInPolygon[1])->coordinates());
        for(int i=1; i<_shortestPathInPolygon.size()-2; i++)
        {
            painter->drawLine(_naivePolygon.vertex(_shortestPathInPolygon[i])->coordinates(),
                              _naivePolygon.vertex(_shortestPathInPolygon[i+1])->coordinates());
        }
        painter->drawLine(_naivePolygon.vertex(_shortestPathInPolygon.end()[-2])->coordinates(),_q->coordinates());
    }
    else
    {
        for (auto elem: _visibilityGraph)
        {
            for(auto vertex: elem.second)
            {
                painter->drawLine(_naivePolygon.vertex(elem.first)->coordinates(),
                                  _naivePolygon.vertex(vertex)->coordinates());
            }
        }
    }
    painter->setPen(regular);
    for (auto v: _naivePolygon.vertices())
    {
        if(!_shortestPathInPolygon.empty())
        {
            int id1 = v->getId();
            int id2 = v->incidentEdge()->twin()->origin()->getId();
            std::vector<int> ids1{id1,id2};
            std::vector<int> ids2{id2,id1};

            auto ind1 = std::search(_shortestPathInPolygon.begin(), _shortestPathInPolygon.end(), ids1.begin(), ids1.end()) != _shortestPathInPolygon.end();
            auto ind2 = std::search(_shortestPathInPolygon.begin(), _shortestPathInPolygon.end(), ids2.begin(), ids2.end()) != _shortestPathInPolygon.end();

            if(ind1 || ind2)continue;
        }
        painter->drawLine(v->coordinates(), v->incidentEdge()->twin()->origin()->coordinates());
    }
}

std::vector<QPointF> SPP::ucitajNasumicneTacke(int brojTacaka) const
{
    std::vector<QPointF> tacke = generisiNasumicneTacke(brojTacaka);

    QPointF maxTacka = tacke[0];

    for (auto i = 1ul; i < tacke.size(); i++) {
        if (tacke[i].x() > maxTacka.x() ||
           (pomocneFunkcije::bliski(tacke[i].x(), maxTacka.x()) && tacke[i].y() < maxTacka.y()))
            maxTacka = tacke[i];
    }

    std::sort(tacke.begin(), tacke.end(), [&](const auto& lhs, const auto& rhs) {
        return pomocneFunkcije::konveksan(maxTacka, lhs, rhs);
    });

    return tacke;
}

std::vector<QPointF> SPP::generisiNasumicneTacke(int brojTacaka) const
{
    static int constexpr DRAWING_BORDER = 10;

    srand(static_cast<unsigned>(time(nullptr)));
    int xMax;
    int yMax;

    if (_pCrtanje)
    {
        xMax = _pCrtanje->width() - DRAWING_BORDER;
        yMax = _pCrtanje->height() - DRAWING_BORDER;
    }
    else
    {
        xMax = CANVAS_WIDTH;
        yMax = CANVAS_HEIGHT;
    }

    int xMin = DRAWING_BORDER;
    int yMin = DRAWING_BORDER;

    std::vector<QPointF> randomPoints;

    int xDiff = xMax-xMin;
    int yDiff = yMax-yMin;
    for(int i=0; i < brojTacaka; i++)
        randomPoints.emplace_back(xMin + rand()%xDiff, yMin + rand()%yDiff);

    return randomPoints;
}

bool SPP::checkDiagonal(Vertex *v, Vertex *u){

    auto v_next = v->incidentEdge()->next()->origin();
    auto v_prev = v->incidentEdge()->prev()->origin();

    if(pomocneFunkcije::konveksan(v_prev->coordinates(), v->coordinates(),v_next->coordinates())){
        if(!pomocneFunkcije::konveksan(v->coordinates(), u->coordinates(), v_prev->coordinates()) ||
           !pomocneFunkcije::konveksan(v->coordinates(), v_next->coordinates(), u->coordinates()))
                 return true;
    } else {
        if(!pomocneFunkcije::konveksan(v->coordinates(), u->coordinates(), v_prev->coordinates()) &&
           !pomocneFunkcije::konveksan(v->coordinates(), v_next->coordinates(), u->coordinates()))
                  return true;
    }

    return false;
}

std::vector<QPointF> SPP::ucitajPodatkeIzDatoteke(std::string imeDatoteke) const
{
    std::ifstream inputFile(imeDatoteke);
    std::vector<QPointF> points;
    int x, y;
    while(inputFile >> x >> y)
        points.emplace_back(x, y);
    return points;
}

Vertex* SPP::generateRandomPositionInPolygon()
{
    int randIndex = (_polygon.vsize()-1)*random()/RAND_MAX;
    return _polygon.vertex(randIndex);
}

bool SPP::seceIvicuPoligona(Vertex* u, Vertex* v)
{
    QPointF presek;
    for(auto k = 0ul; k < _naivePolygon.esize()/2; k++){
        auto edge = _naivePolygon.edge(k);
        if (edge == v->incidentEdge() || edge == v->incidentEdge()->prev() ||
            edge == u->incidentEdge() || edge == u->incidentEdge()->prev())
            continue;

        if (pomocneFunkcije::presekDuzi(QLineF(edge->origin()->coordinates(),
                                               edge->twin()->origin()->coordinates()),
                                        QLineF(v->coordinates(), u->coordinates()),
                                        presek)){
            /* Postoji presek */
            return true;
        }
    }
    return false;
}

void SPP::generateVisibilityGraph()
{
    /* Setting the Id for all vertices in _naivePolygon*/
    for(auto i = 0; i < _naivePolygon.vsize(); i++){
        _naivePolygon.vertex(i)->setId(i);

        //Adding vertex i's two adjacent vertices to the adj matrix to later use for dijkstra
        int V = _naivePolygon.vsize();
        int v1 = i+1 == V ? 0 : i+1;
        int v2 = i-1 == -1 ? V-1 : i-1;

        _visibilityGraph[i].insert(v1);
        _visibilityGraph[i].insert(v2);
    }

    /* Slozenost algoritma: O(n^3) */
    for(auto i = 0; i < _naivePolygon.vsize(); i++){
        auto v = _naivePolygon.vertex(i);
        for (auto j = 0; j < _naivePolygon.vsize(); j++){
            if (i == j)continue;

            auto u = _naivePolygon.vertex(j);
            if(_visibilityGraph[v->getId()].find(u->getId())!=_visibilityGraph[v->getId()].end())continue;

            QPointF presek;

            /* Provera da li je (v, u) spoljasnja dijagonala */
            bool badDiag = checkDiagonal(v, u);

            /* Ako jeste, duz (v, u) nije odgovarajuca, prelazi se na sledecu duz */
            if (badDiag) continue;

            /* Provera da li dijagonala sece neku od ivica poligona */
            for(auto k = 0ul; k < _naivePolygon.esize()/2; k++){
                auto edge = _naivePolygon.edge(k);
                if (edge == v->incidentEdge() || edge == v->incidentEdge()->prev() ||
                    edge == u->incidentEdge() || edge == u->incidentEdge()->prev())
                    continue;

                if (pomocneFunkcije::presekDuzi(QLineF(edge->origin()->coordinates(),
                                                       edge->twin()->origin()->coordinates()),
                                                QLineF(v->coordinates(), u->coordinates()),
                                                presek)){
                    /* Postoji presek */
                    badDiag = true;
                    break;
                }
            }

            /* Ako sece, duz (v, u) nije odgovarajuca, prelazi se na sledecu duz */
            if (badDiag) continue;

            /* Ako ispunjava sve uslove, dodajemo dijagonale (v, u) i (u, v) u graf vidljivosti*/
            if (!badDiag){
                _visibilityGraph[v->getId()].insert(u->getId());
                _visibilityGraph[u->getId()].insert(v->getId());
                AlgoritamBaza_updateCanvasAndBlock()
            }
        }
    }
    AlgoritamBaza_updateCanvasAndBlock()

    for (auto vertexAdj: _visibilityGraph)
    {
        int u = vertexAdj.first;
        for (auto v: _visibilityGraph[u])
        {
            if (_adjMatrix[u][v]!=MAXFLOAT)continue;

            double dist = vDistance(u,v);
            _adjMatrix[u][v] = dist;
            _adjMatrix[v][u] = dist;
        }
    }

    AlgoritamBaza_updateCanvasAndBlock()
    emit animacijaZavrsila();
}

double SPP::vDistance(int u, int v)
{
    auto x1 = _naivePolygon.vertex(u)->coordinates().x();
    auto x2 = _naivePolygon.vertex(v)->coordinates().x();
    auto y1 = _naivePolygon.vertex(u)->coordinates().y();
    auto y2 = _naivePolygon.vertex(v)->coordinates().y();
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
}


