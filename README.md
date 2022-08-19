### Nalazenje puta unutar poligona | Geometrijski algoritmi

Ovaj rad obradjuje metodu odredjivanja najkraćeg puta od jedne poligona do druge, tako da taj put leži u poligonu. 



U radu su implementirane dve verzije algoritma:

1. Naivni
2. Optimizovani



#### Naivni algoritam

Naivni algoritam se sastoji iz tri koraka:

1. Generisanje grafa vidljivosti (`_visibilityGraph`)
2. Pomocu grafa vidljivosti popunjavanje mape susednosti, koja obelezava udaljenost svaka dva cvora u slucaju da su susedna (`_adjMatrix`)
   - u slucaju da se npr cvor `i` i cvor `j` ne vide, u tom slucaju `_adjMatrix[i][j]=MAXFLOAT`
   - u suprotnom, ako se medjusobno vide, tu postavljamo njihovu medjusobnu distancu `_adjMatrix[i][j] = distance(i,j)` (u granicnom slucaju `0` u slucaju da je u pitanju razdaljina cvora od samog sebe)
3. Koriscenje Dijsktrinog algoritma u tandemu sa `_adjMatrix` da bi se utvrdio najkraci put od `_pNaive` do `_qNaive`

Provera vidljivosti dva cvora se proverava tako sto se: provere sledeci uslovi

1. Da li je njihova medjusobna dijagonala unutrasnja ili spoljasnja?
2. Da li linija koja njih povezuje sece bilo koju od stranica poligona?

Najveci nedostatak ovakvog pristupa je slepo izracunavanje vidljivosti izmedju svaka dva cvora, iako nam ona za vecinu cvorova nije relevantna, sto je tim gore zato sto sama operacija izracunavanja vidljivosti izuzetno skupa `O(E)`.



#### Optimizovani algoritam

Optimizovani algoritam je zasnovan na nekoliko koncepata:

1. Pathfinding algoritmu za obilazenje prepreka pomocu grafa vidljivosti
2. Obilazenju lavirinta pravilom leve ruke
3. Path-planning based on visibility graph

Glavni koraci implementiranog algoritma su:

1. Generisanje dve alternativne sekvence koje povezuju krajnje tacke `P` i `Q` preko lanca temena poligona
   - jedan lanac ide u smeru kazaljke na satu, drugi ide suprotno od smera kazaljke na satu
   - ovim delimo cvorove poligona na dve grupe (levu/desnu npr), u zavisnosti od varijaciju problema moze da variraju oba skupa za sve cvorove poligona koji se nalaze u vidokrugu pocetne i krajnje tacke `P` i `Q` (ovo je i neophodna osobina da se resava varijacija problema gde su krajnje tacke lebdece u unutrasnjosti poligona)
2. Svaka od generisanih sekvenci se posmatra kao posebno resenje koje treba "zategnuti"
   - pod zategnuti, impliciram da se u sekvencu uvede najveci moguci broj precica izmedju svaka dva cvora
   - ovo se postize sa dva indeksa `i` i `j`, gde `i` se i krece od `0` do broja cvorova u poligonu `V`, dok se `j` krece u istom rasponu, ali u *suprotnom smeru*
   - ovim poretkom indeksa postizemo da nalazimo najbolje/optimalne precice cim to bude bilo moguce, i necemo pretrazivati nijednu precicu koja bi kasnije bila zamenjena drugom

3. Izracunavamo duzinu obe putanje koje smo dobili i biramo onu kracu
   - s obzirom da je put izmedju `P` i `Q` uvek jedinstven, samo ce jedna od ove dve alternative resenje
   - prosto receno, ne moze doci do toga da optimalni put alterira izmedju dva lanca jer bi se samim tim povecala njegova duzina, jer svaka tacka koju obuhvata levi lanac uvek ima najkracu putanju u "levoj" strani poligona i obrnuto. Njihovo mesanje bi uvek povecalo putanju.  

Precice u ovom algoritmu u realnosti predstavljaju izracunavanje vidljivost izmedju dva cvora. Sto znaci da smo ovim algoritmom dobili izuzetno usmerenu pretragu  koja primenjuje operaciju `isVIsible()` iskljucivo kada je ona neophodna, odnosno kada bi nam doprinela najveci dobitak.