#include "ga07_triangulation.h"
#include <fstream>

Triangulation::Triangulation(QWidget *pCrtanje,
                         int pauzaKoraka,
                         const bool &naivni,
                         std::string imeDatoteke,
                         int brojTacaka)
   : AlgoritamBaza(pCrtanje, pauzaKoraka, naivni),
     _brisucaPravaY(_pCrtanje ? _pCrtanje->height()-3 : 0),
     _statusQueue(StatusQueueCompTriangulation(&_brisucaPravaY))
{
    std::vector<QPointF> tacke;

    /* Ucitavanje tacka */
    if (imeDatoteke != "")
        tacke = ucitajPodatkeIzDatoteke(imeDatoteke);
    else
        tacke = ucitajNasumicneTacke(brojTacaka);

    _monotone = true;
    _polygon.loadData(tacke);
    _naivePolygon.loadData(tacke);
}


void Triangulation::crtajAlgoritam(QPainter *painter) const
{
    if (!painter) return;

    /* Podesavanje razlicitih olovki za crtanje */
    QPen regular = painter->pen();
    regular.setWidth(10);

    QPen magneta = painter->pen();
    magneta.setColor(Qt::magenta);
    magneta.setWidth(2);

    QPen blue = painter->pen();
    blue.setColor(Qt::blue);
    magneta.setWidth(2);

    if (_monotone) {
        /* Crtanje pravougaonika */
        for (auto i=0ul; i < _polygon.vsize() - 1; i++)
            painter->drawLine(_polygon.coordinates(i),
                              _polygon.coordinates(i+1));
        painter->drawLine(_polygon.coordinates(_polygon.vsize()-1),
                          _polygon.coordinates(0));


        /* Crtanje brisuce prave */
        painter->setPen(blue);
        painter->drawLine(QPointF(0, _brisucaPravaY),
                          QPointF(_pCrtanje->width(), _brisucaPravaY));


        /* Crtanje tacaka u redu dogadaja i njihovih tipova. */
        painter->save();
        painter->resetTransform();

        painter->setPen(regular);
        painter->drawText(600, 500, "regular");

        QPen split = regular;
        split.setColor(Qt::darkRed);
        painter->setPen(split);
        painter->drawText(600, 450, "split");

        QPen merge = regular;
        merge.setColor(Qt::darkCyan);
        painter->setPen(merge);
        painter->drawText(600, 400, "merge");

        QPen start = regular;
        start.setColor(Qt::darkYellow);
        painter->setPen(start);
        painter->drawText(600, 350, "start");

        QPen end = regular;
        end.setColor(Qt::darkGreen);
        painter->setPen(end);
        painter->drawText(600, 300, "end");

        painter->restore();

        for(auto v : _eventQueue) {
            switch(v->type())
            {
                case VertexType::START: painter->setPen(start); break;
                case VertexType::REGULAR: painter->setPen(regular); break;
                case VertexType::SPLIT: painter->setPen(split); break;
                case VertexType::END: painter->setPen(end); break;
                case VertexType::MERGE: painter->setPen(merge); break;
            }

            painter->drawPoint(v->coordinates());
        }

        /* Crtanje dijagonala. */
        painter->setPen(magneta);
        for(auto d : _allDiagonals)
            painter->drawLine(d.first->coordinates(), d.second->coordinates());
    } else {
        regular.setWidth(3);
        painter->setPen(regular);
        for (auto f : _polygon.faces()) {
            auto v_start = f->outerComponent();
            auto v_next = v_start;

            if (v_start == nullptr) continue;

            do {
                painter->drawLine(v_next->origin()->coordinates(), v_next->twin()->origin()->coordinates());
                v_next = v_next->next();
            } while (v_next != v_start);
        }

        /* Crtanje brisuce prave */
        painter->setPen(blue);
        painter->drawLine(QPointF(0, _brisucaPravaY),
                          QPointF(_pCrtanje->width(), _brisucaPravaY));

        /* Crtanje dijagonala */
        painter->setPen(magneta);
        for(auto d : _allDiagonals)
            painter->drawLine(d.first->coordinates(), d.second->coordinates());
    }
}

void Triangulation::pokreniAlgoritam()
{
    /* Slozenost algoritma: O(nlogn) */
    initialiseEventQueue();
    monotonePartition();
    _monotone = false;
    connectDiagonalsDCEL();
    _allDiagonals.clear();
    AlgoritamBaza_updateCanvasAndBlock()
    for (auto f : _polygon.faces()) {
        if (f->outerComponent() == nullptr) continue;
        triangulacija(f);
    }
    emit animacijaZavrsila();
}

/**********************************************************************************/
/*                                  Y-MONOTON                                     */
/**********************************************************************************/

void Triangulation::initialiseEventQueue()
{
    for(unsigned i=0; i<_polygon.vsize(); i++){

        auto v = _polygon.vertex(i);
        auto v_sledeci = v->incidentEdge()->twin()->origin();
        auto v_prethodni = v->incidentEdge()->prev()->origin();

    if(pomocneFunkcije::ispod(v_sledeci->coordinates(), v->coordinates()) &&
       pomocneFunkcije::ispod(v_prethodni->coordinates(), v->coordinates())){

        if(pomocneFunkcije::konveksan(v->coordinates(), v_sledeci->coordinates(), v_prethodni->coordinates()))
            v->setType(VertexType::START);
        else v->setType(VertexType::SPLIT);
    }else if(!pomocneFunkcije::ispod(v_sledeci->coordinates(), v->coordinates()) &&
             !pomocneFunkcije::ispod(v_prethodni->coordinates(), v->coordinates())){

        if(pomocneFunkcije::konveksan(v->coordinates(), v_sledeci->coordinates(), v_prethodni->coordinates()))
            v->setType(VertexType::END);
        else v->setType(VertexType::MERGE);

    }else v->setType(VertexType::REGULAR);

        _eventQueue.insert(_polygon.vertex(i));
        AlgoritamBaza_updateCanvasAndBlock()
    }
}

void Triangulation::monotonePartition()
{
    while(!_eventQueue.empty()) {

        auto cvor = *_eventQueue.begin();
        _eventQueue.erase(_eventQueue.begin());

        _brisucaPravaY = cvor->y();
        AlgoritamBaza_updateCanvasAndBlock()

        if(cvor->type() == VertexType::START) {
            handleStartVertex(cvor);
        }
        else if(cvor->type() == VertexType::END) {
            handleEndVertex(cvor);
        }
        else if(cvor->type() == VertexType::MERGE) {
            handleMergeVertex(cvor);
        }
        else if(cvor->type() == VertexType::SPLIT) {
            handleSplitVertex(cvor);
        }
        else {
            handleRegularVertex(cvor);
        }
    }
}

HalfEdge *Triangulation::levaStranica(Vertex *v)
{
    HalfEdge e(v), e_twin(v, &e);
    e.setTwin(&e_twin);

    const auto e_left_it = _statusQueue.lower_bound(&e);
    return e_left_it == _statusQueue.end() ? nullptr : *e_left_it;
}

void Triangulation::handleStartVertex(Vertex *v)
{
    _helpers[v->incidentEdge()] = v;
    _statusQueue.insert(v->incidentEdge());
}

void Triangulation::handleEndVertex(Vertex *v)
{
    if(_helpers[v->incidentEdge()->prev()]->type() == VertexType::MERGE) {
        _allDiagonals.emplace_back(v, _helpers[v->incidentEdge()->prev()]);
        AlgoritamBaza_updateCanvasAndBlock()
    }

    _statusQueue.erase(v->incidentEdge()->prev());
}

void Triangulation::handleSplitVertex(Vertex *v)
{
    /* Pronadi u T stranicu e_left neposredno levo od v. */
    const auto e_left = levaStranica(v);
    if (!e_left) return;

    /* Dodaj dijagonalu v − helper(e_left) u D. */
    _allDiagonals.emplace_back(v, _helpers[e_left]);
    AlgoritamBaza_updateCanvasAndBlock()

    /* Postavi helper(e_left) na v. */
    _helpers[e_left] = v;

    /* Dodaj u T (desnu) stranicu čiji je gornje teme v i postavi njen helper na v. */
    _statusQueue.emplace(v->incidentEdge());
    _helpers[v->incidentEdge()] = v;
}

void Triangulation::handleMergeVertex(Vertex *v)
{
    /* Ako je e stranica čije je donje teme v a poligon joj je lokalno desno,
     * ako je helper(e) merge teme, onda dodaj dijagonalu v − helper(e). */
    auto e = v->incidentEdge()->prev();
    if (_helpers[e]->type() == VertexType::MERGE) {
        _allDiagonals.emplace_back(v, _helpers[e]);
    }

    /* Izbaci e iz T */
    _statusQueue.erase(e);

    /* Pronadi u T stranicu e_left neposredno levo od v */
    const auto e_left = levaStranica(v);
    if (!e_left) return;

    /* Ako je helper(e_left) merge teme, onda dodaj dijagonalu v − helper(e_left). */
    if (_helpers[e_left]->type() == VertexType::MERGE) {
        _allDiagonals.emplace_back(v, _helpers[e_left]);
        AlgoritamBaza_updateCanvasAndBlock();
    }

    /* Postavi helper(e_left) na v. */
    _helpers[e_left] = v;
}

void Triangulation::handleRegularVertex(Vertex *v)
{
    auto v_next = v->incidentEdge()->next()->origin();
    if (pomocneFunkcije::ispod(v_next->coordinates(), v->coordinates())){
        auto helper = _helpers[v->incidentEdge()->prev()];
        if (helper->type()==VertexType::MERGE){
            _allDiagonals.emplace_back(v, helper);
            AlgoritamBaza_updateCanvasAndBlock();
        }
        _statusQueue.erase(v->incidentEdge()->prev());
        _statusQueue.insert(v->incidentEdge());
        _helpers[v->incidentEdge()]=v;
    } else {
        auto left = levaStranica(v);
        if (_helpers[left]->type()==VertexType::MERGE){
            _allDiagonals.emplace_back(v, _helpers[left]);
            AlgoritamBaza_updateCanvasAndBlock();
        }
        _helpers[left]=v;
    }
}

/**********************************************************************************/
/*                             TRIANGULACIJA                                      */
/**********************************************************************************/

void Triangulation::triangulacija(Face *f)
{
    _brisucaPravaY = _pCrtanje->height();
    AlgoritamBaza_updateCanvasAndBlock();

    auto eo = f->outerComponent();
    auto i = eo;

    do{
        _eventQueueTriangulation.insert(i);
        i=i->next();
    }while(i!=eo);

    if (_eventQueueTriangulation.empty()){
        return;
    }
    _stekTriangulacije.push_back(*_eventQueueTriangulation.begin());
    _eventQueueTriangulation.erase(_eventQueueTriangulation.begin());

    if (_eventQueueTriangulation.empty()){
        return;
    }
    _stekTriangulacije.push_back(*_eventQueueTriangulation.begin());
    _eventQueueTriangulation.erase(_eventQueueTriangulation.begin());

    auto u_n= *_eventQueueTriangulation.rbegin();
    for (auto & uj:_eventQueueTriangulation){
    if (uj == *_eventQueueTriangulation.rbegin())
        break;

    if (!istiLanac(_stekTriangulacije.back(), uj)){
        auto uj_1 = _stekTriangulacije.back();
        while(!_stekTriangulacije.empty()){
            auto poslednji = _stekTriangulacije.back();
            _stekTriangulacije.pop_back();

            if (!_stekTriangulacije.empty()){
                _allDiagonals.emplace_back(uj->origin(), poslednji->origin());
                AlgoritamBaza_updateCanvasAndBlock();
            }

            _stekTriangulacije.push_back(uj_1);
            _stekTriangulacije.push_back(uj);
        }

    } else {
          auto poslednji = _stekTriangulacije.back();
          _stekTriangulacije.pop_back();

          while(!_stekTriangulacije.empty()
                &&
                pomocneFunkcije::konveksan(uj->origin()->coordinates(),
                                           poslednji->origin()->coordinates(),
                                           _stekTriangulacije.back()->origin()->coordinates())){
              poslednji = _stekTriangulacije.back();
              _stekTriangulacije.pop_back();
          }

          _stekTriangulacije.push_back(poslednji);
    }
    }

    _stekTriangulacije.pop_back();
    while (!_stekTriangulacije.empty()){
        if (_stekTriangulacije.size()!=1){
        _allDiagonals.emplace_back(u_n->origin(), _stekTriangulacije.back()->origin());
        AlgoritamBaza_updateCanvasAndBlock();
        _stekTriangulacije.pop_back();
        } else{
            _stekTriangulacije.pop_back();
        }
    }
}

bool Triangulation::istiLanac(HalfEdge* e1, HalfEdge* e2)
{
    return (leviLanac(e1, e2) || desniLanac(e1, e2));
}

bool Triangulation::leviLanac(HalfEdge *e1, HalfEdge *e2)
{
    QPointF a = e1->origin()->coordinates();
    QPointF a_next = e1->twin()->origin()->coordinates();

    QPointF b = e2->origin()->coordinates();
    QPointF b_next = e2->twin()->origin()->coordinates();

    return pomocneFunkcije::ispod(a_next, a) &&
           pomocneFunkcije::ispod(b_next, b);
}

bool Triangulation::desniLanac(HalfEdge *e1, HalfEdge *e2)
{
    QPointF a = e1->origin()->coordinates();
    QPointF a_next = e1->twin()->origin()->coordinates();

    QPointF b = e2->origin()->coordinates();
    QPointF b_next = e2->twin()->origin()->coordinates();

    return pomocneFunkcije::ispod(a, a_next) &&
           pomocneFunkcije::ispod(b, b_next);
}


/**********************************************************************************/
/*                          POMOCNE FUNKCIJE                                      */
/**********************************************************************************/

void Triangulation::connectDiagonalsDCEL()
{
    std::map<Vertex *, std::set<HalfEdge *, DiagonalsAddDECELComp>> allDiagonals;
    std::map<Vertex *, std::set<HalfEdge *, DiagonalsAddDECELComp>>::iterator it;


    for(auto &pair : _allDiagonals) {
        auto ei1 = new HalfEdge(pair.first);
        auto ei2 = new HalfEdge(pair.second, ei1);
        ei1->setTwin(ei2);

        _polygon.insertEdge(ei1);
        _polygon.insertEdge(ei2);

        it = allDiagonals.find(pair.first);
        if (it != allDiagonals.end()) {
            allDiagonals[pair.first].emplace(ei1);
        } else {
            std::set<HalfEdge *, DiagonalsAddDECELComp> s1(DiagonalsAddDECELComp(pair.first));
            s1.emplace(ei1);
            allDiagonals[pair.first] = s1;
        }

        it = allDiagonals.find(pair.second);
        if (it != allDiagonals.end()) {
            allDiagonals[pair.second].emplace(ei2);
        } else {
            std::set<HalfEdge *, DiagonalsAddDECELComp> s2(DiagonalsAddDECELComp(pair.second));
            s2.emplace(ei2);
            allDiagonals[pair.second] = s2;
        }
    }

    auto f_old = _polygon.face(1);
    for(auto &v : allDiagonals) {
        auto v_prev = v.first->incidentEdge()->prev();
        auto v_next = v.first->incidentEdge();

        for(auto &d : v.second) {
            v_prev->setNext(d);
            d->setPrev(v_prev);

            v_prev = d->twin();

            /* Dodavanje polja (face) ne radi kako treba.
             * Ponekad pravi greske i za outerComponent za face se zadaje pogresne
             * ivice sto dovodi do greske. */
            if (d->incidentFace() == nullptr) {
                auto f = new Face();
                _polygon.insertFace(f);

                if (sameDirectionVectors(f_old->outerComponent(), d)) {
                    f->setOuterComponent(d);
                    f_old->setOuterComponent(d->twin());
                    d->setIncidentFace(f);
                    d->twin()->setIncidentFace(f_old);
                } else {
                    f->setOuterComponent(d->twin());
                    f_old->setOuterComponent(d);
                    d->twin()->setIncidentFace(f);
                    d->setIncidentFace(f_old);
                }

                f_old = f;
            }
        }
        v_prev->setNext(v_next);
        v_next->setPrev(v_prev);
    }
}

bool Triangulation::sameDirectionVectors(HalfEdge *e1, HalfEdge *e2)
{
    auto a = e1->origin();
    auto b = e1->twin()->origin();

    auto c = e2->origin();
    auto d = e2->twin()->origin();

    return  (b->x() - a->x()) *
            (d->x() - c->x()) +
            (b->y() - a->y()) *
            (d->y() - c->y()) > 0;
}

std::vector<QPointF> Triangulation::ucitajPodatkeIzDatoteke(std::string imeDatoteke) const
{
    std::ifstream inputFile(imeDatoteke);
    std::vector<QPointF> points;
    int x, y;
    while(inputFile >> x >> y)
        points.emplace_back(x, y);
    return points;
}

std::vector<QPointF> Triangulation::generisiNasumicneTacke(int brojTacaka) const
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

std::vector<QPointF> Triangulation::ucitajNasumicneTacke(int brojTacaka) const
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

/**********************************************************************************/
/*                         NAIVNI ALGORITAM                                       */
/**********************************************************************************/

void Triangulation::pokreniNaivniAlgoritam()
{
    /* Slozenost algoritma: O(n^3) */
    for(auto i = 0ul; i < _naivePolygon.vsize(); i++){
        auto v = _naivePolygon.vertex(i);

        for (auto j = i+1; j < _naivePolygon.vsize(); j++){
            auto u = _naivePolygon.vertex(j);

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

            /* Provera da li (v, u) sece neku od trenutnih dijagonala */
            for (auto diag: _naiveDiagonals){
                    if (diag.first == v || diag.second == v ||
                        diag.first == u || diag.second == u)
                        continue;

                    if (pomocneFunkcije::presekDuzi(QLineF(diag.first->coordinates(),
                                                           diag.second->coordinates()),
                                                    QLineF(v->coordinates(),u->coordinates()),
                                                    presek)){
                        /* Sece dijagonalu */
                        badDiag = true;
                        break;
                    }
             }

            /* Ako ispunjava sve uslove, dodajemo dijagonalu (v, u) */
            if (!badDiag){
                _naiveDiagonals.emplace_back(v,u);
                AlgoritamBaza_updateCanvasAndBlock()
            }
        }
    }

    AlgoritamBaza_updateCanvasAndBlock()
    emit animacijaZavrsila();
}

bool Triangulation::checkDiagonal(Vertex *v, Vertex *u){

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

void Triangulation::crtajNaivniAlgoritam(QPainter *painter) const
{
    if (!painter) return;

    QPen regular = painter->pen();
    regular.setWidth(3);

    QPen green = painter->pen();
    green.setColor(Qt::darkGreen);
    green.setWidth(2);
    painter->setPen(green);

    for(auto diag: _naiveDiagonals)
        painter->drawLine(diag.first->coordinates(),diag.second->coordinates());

    painter->setPen(regular);
    for (auto v: _naivePolygon.vertices())
        painter->drawLine(v->coordinates(), v->incidentEdge()->twin()->origin()->coordinates());
}

const DCEL &Triangulation::getNaivePolygon() const
{
    return _naivePolygon;
}

const std::vector<std::pair<Vertex *, Vertex *> > &Triangulation::getNaiveDiagonals() const
{
    return _naiveDiagonals;
}
