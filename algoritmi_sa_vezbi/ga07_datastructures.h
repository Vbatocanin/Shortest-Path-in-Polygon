#ifndef GA07_DATASTRUCTURES_H
#define GA07_DATASTRUCTURES_H

#include "ga06_dcel.h"
#include "pomocnefunkcije.h"

///
/// \brief The EventQueueComp struct
/// Brisuca prava se krece od jedne ka drugoj EventTacki
///
struct EventQueueCompTriangulation
{
    bool operator()(const Vertex* levi, const Vertex* desni) const
    {
        return levi->coordinates().y() > desni->coordinates().y() ||
              (pomocneFunkcije::bliski(levi->coordinates().y(), desni->coordinates().y()) &&
                                       levi->coordinates().x() < desni->coordinates().y());
    }
};

///
/// \brief The StatusQueueCompTriangulation struct
///
/// Status za Monotono Particionisanje
///
struct StatusQueueCompTriangulation
{
private:
    const double *const yBrisucaPrava;

public:
    StatusQueueCompTriangulation(double *y)
        : yBrisucaPrava(y)
    {}

    bool operator()(const HalfEdge* line1, const HalfEdge* line2) const
    {

       QLineF brisucaPrava(0, *yBrisucaPrava, 5, *yBrisucaPrava);
       QPointF presek1;
       QPointF presek2;

       QLineF line11(line1->origin()->coordinates(), line1->twin()->origin()->coordinates());
       QLineF line22(line2->origin()->coordinates(), line2->twin()->origin()->coordinates());

       pomocneFunkcije::presekDuzi(brisucaPrava, line11, presek1);
       pomocneFunkcije::presekDuzi(brisucaPrava, line22, presek2);

       return presek1.x() < presek2.x();
    }
};


struct DiagonalsAddDECELComp
{
private:
    Vertex* _vertex;

public:
    DiagonalsAddDECELComp(Vertex* vertex = nullptr)
        :_vertex(vertex)
    {}

    bool operator()(const HalfEdge* line1, const HalfEdge* line2) const
    {
        Vertex* v1 = line1->origin() != _vertex ? line1->origin() :
                                                  line1->twin()->origin();

        Vertex* v2 = line2->origin() != _vertex ? line2->origin() :
                                                  line2->twin()->origin();

        return !pomocneFunkcije::konveksan(_vertex->coordinates(), v1->coordinates(), v2->coordinates());
    }
};

///
/// \brief The EventQueueCompTriangulation2 struct
/// U Event za Triangulaciju pamtimo HalfEdge, a ne vertex
/// (kao sto je dato u knjizi)
/// zbog problema u implementaciji koji postoje ako se koristi vertex
///
/// TODO
struct EventQueueCompTriangulation2
{
    bool operator()(const HalfEdge* lhs, const HalfEdge* rhs) const
    {
        return lhs->origin()->y() > rhs->origin()->y() ||
              (pomocneFunkcije::bliski(lhs->origin()->y(), rhs->origin()->y()) &&
                                       lhs->origin()->x() < rhs->origin()->x());

    }
};


#endif // GA07_DATASTRUCTURES_H
