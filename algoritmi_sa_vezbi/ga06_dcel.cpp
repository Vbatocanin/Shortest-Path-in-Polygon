#include "ga06_dcel.h"

#include <iostream>
#include <fstream>
#include "config.h"
#include <cstdlib>
#include <iostream>
#include <math.h>
bool circumCircleContains(Vertex *a, Vertex* b, Vertex* c,Vertex* v);
bool isIllegal(Vertex* v, HalfEdge *e);
Vertex* returnPK(Vertex* v, HalfEdge* e);
Face* Face::CopyOfTriangle(Face* triangle){
    HalfEdge* e1_old = triangle->outerComponent();
    HalfEdge* e2_old = e1_old->next();
    HalfEdge* e3_old = e2_old->next();

    Vertex* a = e1_old->origin();
    Vertex* b = e2_old->origin();
    Vertex* c = e3_old->origin();

    HalfEdge* e1 = new HalfEdge(); HalfEdge* e1_twin = new HalfEdge();
    HalfEdge* e2 = new HalfEdge(); HalfEdge* e2_twin = new HalfEdge();
    HalfEdge* e3 = new HalfEdge(); HalfEdge* e3_twin = new HalfEdge();

    /*Vertex* v1_new = new Vertex(a->x, a->y, e1);
    Vertex* v2_new = new Vertex(b->x, b->y, e2);
    Vertex* v3_new = new Vertex(c->x, c->y, e3);
    v1_new->id = a->id;
    v2_new->id = b->id;
    v3_new->id = c->id;

    Face* d = new Face(e1);*/
  /*  e1->setEdgeParameters(v1_new, e1_twin, d, e2, e3);
    e2->setEdgeParameters(v2_new, e2_twin, d, e3, e1);
    e3->setEdgeParameters(v3_new, e3_twin, d, e1, e2);

    e1_twin->setEdgeParameters(v2_new, e1, d, e3_twin, e2_twin);
    e2_twin->setEdgeParameters(v3_new, e2, d, e1_twin, e3_twin);
    e3_twin->setEdgeParameters(v1_new, e3, d, e2_twin, e1_twin);

    return d;*/
}

int min(int x, int y){
    return x<y? x:y;
}

Vertex* returnPK(Vertex* v, HalfEdge* e){
    if (e->isLeft(v)){
        return e->twin()->prev()->origin();
    } else{
        return e->prev()->origin();
    }
}


bool isIllegal(Vertex* v, HalfEdge *e){
    Vertex* p_k = returnPK(v, e);

    if (e->isBordering())
         return false;

    if (p_k->getId()>0 && v->getId()>0 && e->origin()->getId()>0
            && e->twin()->origin()->getId()>0)
        return circumCircleContains(p_k,e->origin(), e->twin()->origin(),
                                    v);
    if (min(p_k->getId(),v->getId())>=min(e->origin()->getId(),
                                          e->twin()->origin()->getId()))
        return true;
}

void DCEL::legalize_edge(Vertex* v, HalfEdge* e){
    if (isIllegal(v,e)){
        Vertex* p_k = returnPK(v,e);
         if (e->isLeft(v)){
            legalize_edge(v, e->twin()->next());
            legalize_edge(v, e->twin()->prev());
         } else{

            legalize_edge(v, e->next());
            legalize_edge(v, e->prev());
        }
        flip(v, e, p_k);
    }
}

void DCEL::flip(Vertex* v, HalfEdge* e, Vertex* p_k){
    if (e->isLeft(v)){
        HalfEdge * nova_ivica = new HalfEdge(v);
        e->next()->setNext(nova_ivica);
        nova_ivica->setPrev(e->next());
        // cd je previous od da
        e->next()->setPrev(e->twin()->prev());
        nova_ivica->setNext(e->twin()->prev());
        e->twin()->prev()->setPrev(nova_ivica);
        HalfEdge * nova_ivica_twin = new HalfEdge(p_k, nova_ivica);
        // setovati i za twin odg.
    } else {

    }
}
double Vertex::dist2(const Vertex& v) const{
    return (x()-v.x())*(x()-v.x()) + (y()-v.y())*(y()-v.y());
}
double Vertex::norm2() const{
    return this->x()*this->x() + this->y()*this->y();
}

bool Vertex::operator<(const Vertex &other) const
{
    return this->x() < other.x() || (this->x() == other.x() && this->y() < other.y());
}
bool circumCircleContains(Vertex *a, Vertex* b, Vertex* c,Vertex* v){
    const double ab = a->norm2();
    const double cd = b->norm2();
    const double ef = c->norm2();

    const double ax = a->x();
    const double ay = a->y();
    const double bx = b->x();
    const double by = b->y();
    const double cx = c->x();
    const double cy = c->y();

    const double circum_x = (ab*(cy-by)+cd*(ay-cy)+ef*(by-ay))/(ax*(cy-by)+bx*(ay-cy)+cx*(by-ay));
    const double circum_y = (ab*(cx-bx)+cd*(ax-cx)+ef*(bx-ax))/(ay*(cx-bx)+by*(ax-cx)+cy*(bx-ax));
    Vertex* circum = new Vertex(circum_x, circum_y);
    const double circum_radius = a->dist2(*circum);
    const double dist = v->dist2(*circum);
    return dist <= circum_radius;

}
bool HalfEdge::containsPoint(Vertex *v)const{
    Vertex* a = this->origin();
    Vertex* b = this->twin()->origin();

    if (fabs(b->coordinates().x()-a->coordinates().x())< EPS){
        if (fabs(b->coordinates().y()-a->coordinates().y())<EPS)
            return ((fabs(b->coordinates().x()-v->coordinates().x())<EPS)
                    && (fabs(b->coordinates().y()-v->coordinates().y())<EPS));
        else
            return (fabs(b->coordinates().x()-v->coordinates().x())<EPS)
                    &&
                    (((a->coordinates().y()<v->coordinates().y())
                    && (v->coordinates().y()<b->coordinates().y()))
                    ||  ((b->coordinates().y()<v->coordinates().y())
                         && (v->coordinates().y()<a->coordinates().y())));

    }
    double x_rez =(v->coordinates().x()-a->coordinates().x())/(b->coordinates().x()
                                                               -a->coordinates().x());
    double y_rez =(v->coordinates().y()-a->coordinates().y())/
            (b->coordinates().y()-a->coordinates().y());
/*
    if (fabs(x_rez-y_rez)<epsilon)
        return true;
    return false;*/
}

bool HalfEdge::isLeft(Vertex* v) const{
    Vertex* e_from = this->origin();
    Vertex* e_to = this->twin()->origin();

    if ((e_from->getId()==-2) && (e_to->getId()==-1))
        return true;

    if ((e_from->getId() == -1) || (e_to->getId() == -2)){
        return v < e_to;
    } else if ((e_from->getId() == -2) || (e_to->getId() == -1)){
        return v > e_to;
    } else {
        double cross_prod=(e_to->coordinates().x() - e_from->coordinates().x())*
                (v->coordinates().y() - e_from->coordinates().y())
                -(e_to->coordinates().y() - e_from->coordinates().y())
                *(v->coordinates().x()-e_from->coordinates().x());
        return cross_prod < EPS;
    }
}

bool Face::inTriangle(Vertex* v) const{
    HalfEdge* e1 = this->outerComponent();
    HalfEdge* e2 = e1->next();
    HalfEdge* e3 = e2->next();

    if (e1->isBordering() && e2->isBordering() && e3->isBordering())
        return true;

    return e1->isLeft(v) && e2->isLeft(v) && e3->isLeft(v);
}


bool HalfEdge::isBordering() const{
    Vertex* from = this->origin();
    Vertex* to = this->twin()->origin();

    if ((from->getId() >=-2 && from->getId()<=0) and (to->getId() >=-2 && to->getId()<=0))
        return true;
    return false;
}
void DCEL::add_point_in_triangle(Face * triangle, Vertex *v){
        HalfEdge* e1 = triangle->outerComponent();
        HalfEdge* e2 = e1->next();
        HalfEdge* e3 = e2->next();

        Vertex* a = e1->origin();
        Vertex* b = e2->origin();
        Vertex* c = e3->origin();


        HalfEdge* a1 = new HalfEdge(); HalfEdge* a2 = new HalfEdge();
        HalfEdge* b1 = new HalfEdge(); HalfEdge* b2 = new HalfEdge();
        HalfEdge* c1 = new HalfEdge(); HalfEdge* c2 = new HalfEdge();

        v->setIncidentEdge(a2);

        Face* f1 = new Face(a2, std::vector<HalfEdge*>{});
        Face* f2 = new Face(b2, std::vector<HalfEdge*>{});
        Face* f3 = new Face(c2, std::vector<HalfEdge*>{});

        /*
        e1->next = b1; e1->prev=a2; e1->incidentFace=f1;
        e2->next = c1; e2->prev=b2; e2->incidentFace=f2;
        e3->next = a1; e3->prev=c2; e3->incidentFace=f3;
        */
        a1->setEdgeParameters(a,a2,f3,c2,e3);
        a2->setEdgeParameters(v,a1,f1,e1,b1);
        b1->setEdgeParameters(b,b2,f1,a2,e1);
        b2->setEdgeParameters(v,b1,f2,e2,c1);
        c1->setEdgeParameters(c,c2,f2,b2,e2);
        c2->setEdgeParameters(v,c1,f3,e3,a1);

        (triangle->getChildrenFaces()).push_back(f1);
        (triangle->getChildrenFaces()).push_back(f2);
        (triangle->getChildrenFaces()).push_back(f3);

        /*
        f1->printFace();
        f2->printFace();
        f3->printFace();*/
}
void DCEL::add_point_on_edge(Face * triangle, HalfEdge* e, Vertex *v){
        HalfEdge* e_next = e->next();
        HalfEdge* e_prev = e->prev();
        HalfEdge* e_twin = e->twin();

        HalfEdge* e1 = new HalfEdge();
        HalfEdge* e1_twin = new HalfEdge();
        HalfEdge* e2 = new HalfEdge();
        HalfEdge* e2_twin = new HalfEdge();
        HalfEdge* e3 = new HalfEdge();
        HalfEdge* e3_twin = new HalfEdge();
        HalfEdge* e4 = new HalfEdge();
        HalfEdge* e4_twin = new HalfEdge();

        Face* f1 = new Face(e1,std::vector<HalfEdge*>{});
        Face* f2 = new Face(e4,std::vector<HalfEdge*>{});
        Face* f3 = new Face(e2,std::vector<HalfEdge*>{});
        Face* f4 = new Face(e2_twin,std::vector<HalfEdge*>{});

        v->setIncidentEdge(e2);

        e1->setEdgeParameters(e->origin(), e1_twin, f1, e3, e_prev);
        e1_twin->setEdgeParameters(v, e1, f2, e_twin->next(), e4);
        e2->setEdgeParameters(v, e2_twin, f3, e_next, e3_twin);
        e2_twin->setEdgeParameters(e_twin->origin(), e2, f4, e4_twin, e_twin->prev());
        e3->setEdgeParameters(v, e3_twin, f1, e_prev, e1);
        e3_twin->setEdgeParameters(e_prev->origin(), e3, f3, e2, 	e_next);
        e4->setEdgeParameters((e_twin->prev())->origin(), e4_twin, f2, e1_twin,e_twin->next());
        e4_twin->setEdgeParameters(v, e4_twin, f4, e_twin->prev(), e2_twin);

        /* setovati spoljne ivice, njihove parametre*/

        (triangle->getChildrenFaces()).push_back(f1);
        (triangle->getChildrenFaces()).push_back(f3);
        ((e_twin->incidentFace())->getChildrenFaces()).push_back(f2);
        ((e_twin->incidentFace())->getChildrenFaces()).push_back(f4);

        /*f1->printFace();
        f2->printFace();
        f3->printFace();
        f4->printFace();*/
}

void DCEL::add_point_general(Face* triangle, Vertex* v){
    HalfEdge* e1 = triangle->outerComponent();
    HalfEdge* e2 = e1->next();
    HalfEdge* e3 = e2->next();


    if (e1->containsPoint(v)){
        add_point_on_edge(triangle, e1, v);
    } else if (e2->containsPoint(v)){
        add_point_on_edge(triangle, e2, v);
    } else if (e3->containsPoint(v)){
        add_point_on_edge(triangle, e3, v);
    } else{
        add_point_in_triangle(triangle, v);
    }

}
/***********************************************************************/
/*                               DCEL                                  */
/***********************************************************************/
DCEL::DCEL(std::string imeDatoteke, int h, int w)
{
    std::ifstream in(imeDatoteke);
    std::string tmp;

    in >> tmp;
    if (tmp.compare("OFF") != 0) {
        std::cout << "Wrong file format: " << imeDatoteke << std::endl;
        exit(EXIT_FAILURE);
    }

    unsigned vertexNum, edgeNum, faceNum;
    in >> vertexNum >> faceNum >> edgeNum;

    /* citamo teme po teme iz fajla
     * i pamtimo ga u nizu temena */
    float tmpx, tmpy, tmpz;
    float x, y;
    Vertex* v;
    for(auto i=0ul; i<vertexNum; i++) {
        in >> tmpx >> tmpy >> tmpz;
        x = (tmpx + 1) / 2 * w;
        y = (tmpy + 1) / 2 * h;
        v = new Vertex(x, y);
        _vertices.push_back(v);
    }

    /* za svaki poligon (polje)
     * za svako teme tog poligona napravimo polustranicu sa tim temenom kao pocetkom
     * zatim prodjemo kroz napravljene polustranice i popunimo prethodnu i sledecu */
    for(auto i=0ul; i<faceNum; i++) {
        auto f = new Face();
        unsigned broj_temena;
        in >> broj_temena;
        std::vector<HalfEdge*> edges;
        for(auto j = 0ul; j < broj_temena; ++j){
            unsigned indeksTemena;
            in >> indeksTemena;
            HalfEdge *new_halfedge = new HalfEdge(_vertices[indeksTemena]);
            edges.push_back(new_halfedge);
            _edges.push_back(new_halfedge);
            _vertices[indeksTemena]->setIncidentEdge(new_halfedge);
        }
        for(auto j = 0ul; j < broj_temena; j++){
            edges[j]->setNext(edges[(j + 1) % broj_temena]);
            edges[j]->setPrev(edges[(j - 1 + broj_temena) % broj_temena]);
            edges[j]->setIncidentFace(f);
        }
        f->setOuterComponent(edges[0]);
        _faces.push_back(f);
    }

    /* mapa u kojoj mozemo naci stranicu na osnovu temena na koji "pokazuje"
     * u njoj cemo pamtiti spoljasnje ivice i koristimo je da bi ih kasnije povezali */
    std::map<Vertex*, HalfEdge*> spoljasnje_ivice;

    auto spoljasnost = new Face();

    /* za svaku polustranicu AB pokusavamo da nadjemo polustranicu BA
     * ako takva ne postoji to znaci da smo naisli na spoljasnju polustranicu
     * i nju smestamo u mapu spoljasnje ivice za kasniju obradu */
    for(auto edge : _edges){
        auto uslov = [=](HalfEdge* e){
            return e->origin() == edge->next()->origin()
                   && e->next()->origin() == edge->origin();
        };
        auto twin = std::find_if(std::begin(_edges), std::end(_edges), uslov);
        if(twin == std::end(_edges)){
            HalfEdge *outer_edge = new HalfEdge();
            edge->setTwin(outer_edge);
            outer_edge->setIncidentFace(spoljasnost);
            outer_edge->setTwin(edge);
            outer_edge->setOrigin(edge->next()->origin());
            spoljasnje_ivice[edge->origin()] = outer_edge;
        }
        else{
            edge->setTwin(*twin);
        }
    }

    /* za svaku spoljasnju ivicu trazimo spoljasnju koja joj prethodi
     * to mozemo jednostavno uraditi zbog nacina pravljenja mape spoljasnje_ivice
     * nakon sto zavrsimo "sredjivanje" te poluivice dodajemo je u listu poluivica DCEL strukture */
    for(auto outer_edge : spoljasnje_ivice){
        outer_edge.second->setPrev(spoljasnje_ivice[outer_edge.second->origin()]);
        spoljasnje_ivice[outer_edge.second->origin()]->setNext(outer_edge.second);
        _edges.push_back(outer_edge.second);
    }

    spoljasnost->setInnerComponent(spoljasnje_ivice.begin()->second);
    _faces.push_back(spoljasnost);
}

DCEL::DCEL(const std::vector<QPointF> &tacke)
{
    loadData(tacke);
}

DCEL::~DCEL() {
    for (auto v: _vertices)
        delete v;

    for (auto e: _edges)
        delete e;

    for (auto f: _faces)
        delete f;
}

size_t DCEL::vsize() const
{
    return _vertices.size();
}

Vertex *DCEL::vertex(size_t i) const
{
    return _vertices[i];
}

const QPointF &DCEL::coordinates(size_t i) const
{
    return _vertices[i]->coordinates();
}

const std::vector<Vertex *> &DCEL::vertices() const
{
    return _vertices;
}

size_t DCEL::esize() const
{
    return _edges.size();
}

HalfEdge *DCEL::edge(size_t i) const
{
    return _edges[i];
}


const std::vector<HalfEdge *> &DCEL::edges() const
{
    return _edges;
}

size_t DCEL::fsize() const
{
    return _faces.size();
}

Face *DCEL::face(size_t i) const
{
    return _faces[i];
}


const std::vector<Face *> &DCEL::faces() const
{
    return _faces;
}


void DCEL::setFaces(const std::vector<Face *> &faces)
{
    _faces = faces;
}

void DCEL::loadData(const std::vector<QPointF> &tacke)
{
    /* tacke moraju biti zadate u smeru suprotnom kazaljci na satu
       inace ova metoda nije validna */

    auto inner = new Face();
    auto outer = new Face();

    _faces.push_back(outer);
    _faces.push_back(inner);

    for(auto &tacka : tacke) {
        auto newVertex = new Vertex(tacka);
        _vertices.push_back(newVertex);

        auto newEdge = new HalfEdge(newVertex);
        newEdge->setIncidentFace(inner);
        _edges.push_back(newEdge);
    }

    for(auto j = 0ul; j < _vertices.size(); ++j) {
        /* Postavljanje odgovarajuce ivice za Vertex */
        _vertices[j]->setIncidentEdge(_edges[j]);

        /* Postavljanje next ivice */
        if(j != _vertices.size()-1)
            _edges[j]->setNext(_edges[j+1]);
        else
            _edges[j]->setNext(_edges[0]);

        /* Postavljanje prev ivice */
        if(j != 0)
            _edges[j]->setPrev(_edges[j-1]);
        else
            _edges[j]->setPrev(_edges[_vertices.size()-1]);

        /* Postavljamo twin ivice */
        auto twinEdge = new HalfEdge(_vertices[(j+1) % _vertices.size()],
                                     _edges[j], nullptr, nullptr, outer);
        _edges[j]->setTwin(twinEdge);
        _edges.push_back(twinEdge);
        /* Sada smo poceli da stavljamo twin ivice, one pocinju od pozicije n+j u vektoru */
    }

    for(auto k = _vertices.size(); k < _edges.size(); ++k) {
        /* Postavljanje next ivice za twin */
        if(k != _edges.size()-1)
            _edges[k]->setNext(_edges[k+1]);
        else
            _edges[k]->setNext(_edges[_vertices.size()]);

        /* Postavljanje prev ivice za twin */
        if(k != 0)
            _edges[k]->setPrev(_edges[k-1]);
        else
            _edges[k]->setPrev(_edges[_edges.size()-1]);
    }

    /* Ovde je izabran _edges[0] umesto njegovog twin-a, zato sto smo izabrali
     * da radimo u smeru suprotnom kazaljci na satu */
    inner->setOuterComponent(_edges[0]);

    outer->setInnerComponent(_edges[0]->twin());
    outer->setOuterComponent(nullptr); /* Outer nema granicu */
}

void DCEL::insertEdge(HalfEdge *e)
{
    _edges.push_back(e);
}

void DCEL::insertFace(Face *f)
{
    _faces.push_back(f);
}

HalfEdge *DCEL::findEdge(Vertex *start, Vertex *end)
{
    for(HalfEdge* e : _edges)
    {
        if (e->origin() == start && e->twin()->origin() == end)
            return e;
    }

    return nullptr;
}
/*
***********************************************************************
*                             VERTEX                                  *
***********************************************************************
*/
Vertex::Vertex()
    : _incidentEdge{nullptr}
{}

Vertex::Vertex(float x, float y)
    : _coordinates{static_cast<qreal>(x), static_cast<qreal>(y)}
{}

Vertex::Vertex(const QPointF &coordinates, HalfEdge *incidentEdge)
    : _coordinates{coordinates}, _incidentEdge{incidentEdge}
{}

qreal Vertex::x() const
{
    return _coordinates.x();
}

qreal Vertex::y() const
{
    return _coordinates.y();
}

const QPointF &Vertex::coordinates() const
{
    return _coordinates;
}

void Vertex::setCoordinates(const QPointF &coordinates)
{
    _coordinates = coordinates;
}

HalfEdge *Vertex::incidentEdge() const
{
    return _incidentEdge;
}

void Vertex::setIncidentEdge(HalfEdge *incidentEdge)
{
    _incidentEdge = incidentEdge;
}
/*
***********************************************************************
*                           HALFEDGE                                 *
***********************************************************************
*/

HalfEdge::HalfEdge()
    : _origin{nullptr}, _twin{nullptr}, _next{nullptr}, _prev{nullptr}, _incidentFace{nullptr}
{}

HalfEdge::HalfEdge(Vertex *origin, HalfEdge *twin, HalfEdge *next, HalfEdge *prev, Face *incidentFace)
    : _origin(origin), _twin(twin), _next(next), _prev(prev), _incidentFace(incidentFace)
{}

Vertex *HalfEdge::origin() const
{
    return _origin;
}

void HalfEdge::setOrigin(Vertex *origin)
{
    _origin = origin;
}

HalfEdge *HalfEdge::twin() const
{
    return _twin;
}

void HalfEdge::setTwin(HalfEdge *twin)
{
    _twin = twin;
}

HalfEdge *HalfEdge::next() const
{
    return _next;
}

void HalfEdge::setNext(HalfEdge *next)
{
    _next = next;
}

HalfEdge *HalfEdge::prev() const
{
    return _prev;
}

void HalfEdge::setPrev(HalfEdge *prev)
{
    _prev = prev;
}

Face *HalfEdge::incidentFace() const
{
    return _incidentFace;
}

void HalfEdge::setIncidentFace(Face *incidentFace)
{
    _incidentFace = incidentFace;
}

/*
***********************************************************************
*                             FACE (FIELD)                            *
***********************************************************************
*/
Face::Face()
    :_outerComponent{nullptr}
{}

Face::Face(HalfEdge *outerComponent, const std::vector<HalfEdge *> &innerComponent)
    :_outerComponent{outerComponent}, _innerComponents{innerComponent}
{}

HalfEdge *Face::outerComponent() const
{
    return _outerComponent;
}

void Face::setOuterComponent(HalfEdge *outerComponent)
{
    _outerComponent = outerComponent;
}

const std::vector<HalfEdge *> &Face::innerComponents() const
{
    return _innerComponents;
}

HalfEdge *Face::innerComponent() const
{
    return _innerComponents.empty() ? nullptr : _innerComponents.front();
}

void Face::setInnerComponents(const std::vector<HalfEdge *> &innerComponents)
{
    _innerComponents = innerComponents;
}

void Face::setInnerComponent(HalfEdge* innerComponent)
{
    _innerComponents.push_back(innerComponent);
}
