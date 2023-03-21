#include <vector>
#include <cmath>
#include <iostream>

const double pi = 3.14159265358;
const double acc = 1e-6;

class Line;

struct Point {
    Point() {};
    Point(const double& x, const double& y) : x(x), y(y) {}

    bool operator==(const Point& another) const { return (abs(x - another.x) < acc) && (abs(y - another.y) < acc); }

    bool operator!=(const Point& another) const { return !(*this == another); }

    double dist(const Point& p) const { return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y)); }

    void scale(double a);
    void rotate(double a);
    void reflect(const Point& p);
    void reflect(const Line& l);

    double x = 0;
    double y = 0;
};

class Line {
public:

    Line() : Line(0,0) {}

    Line(double a, double b, double c) : a(a), b(b), c(c) {}
    Line(const Point& a, const Point& b);
    Line(const Point& p, const double& a) : Line(a,p.y - (p.x * a)) {}
    Line(const double& k, const double& b) : Line(k, -1, b) {}



    bool operator==(const Line& another) const;

    bool operator!=(const Line& another) const { return !(*this == another);}

    Line& rotate(const Point& p, const double& d);

    //bool is_par(const Line& l) const { return abs((l.b * a) - (l.a * b)) < acc; }

    Point cross(const Line& l) const;
    Line normal(const Point& p) const;

    double a = 0;
    double b = 0;
    double c = 0;
};

class Shape {
public:
    virtual void rotate(Point center, double angle) = 0;
    virtual void reflect(Point center) = 0;
    virtual void reflect(Line axis) = 0;
    virtual void scale(Point center, double coefficient) = 0;

    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool operator==(const Shape& another) const = 0;
    virtual bool operator!=(const Shape& another) const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual bool containsPoint(Point point) const = 0;
    virtual ~Shape() = 0;
};

class Ellipse : public Shape {
public:

    Ellipse(const Point& f1, const Point& f2, double r) : f1(f1), f2(f2), r(r) {}

    virtual ~Ellipse() override {}

    std::pair<Point,Point> focuses() const {return {f1,f2};}

    std::pair<Line, Line> directrices() const;

    double eccentricity() const { return (c() / a()); }
    Point center() const { return Point((f1.x + f2.x) / 2, (f1.y + f2.y) / 2); }

    virtual void rotate(Point center, double angle) override;
    virtual void reflect(Point center) override;
    virtual void reflect(Line axis) override;
    virtual void scale(Point center, double coefficient) override;

    virtual double perimeter() const override { return 4 * a() * std::comp_ellint_2(eccentricity()); }
    virtual double area() const override{ return pi * a() * b(); }

    virtual bool operator==(const Shape& other) const override;
    virtual bool operator!=(const Shape& other) const override;

    virtual bool isCongruentTo(const Shape& another) const override;
    virtual bool isSimilarTo(const Shape& another) const override;
    virtual bool containsPoint(Point point) const override{ return f1.dist(point) + f2.dist(point) < acc + r; }

protected:
    double a() const { return r / 2; }
    double b() const { return sqrt((a() * a()) - (c() * c())); }
    double c() const { return (f1.dist(f2) / 2); }

    Point f1;
    Point f2;
    double r = 0;
};

class Circle : public Ellipse {
public:
    Circle(const Point& p, double r) : Ellipse(p, p, r * 2) {}

    virtual ~Circle() {};

    double radius() const { return r/2; }

};

class Polygon : public Shape {
public:
    Polygon();
    Polygon(const std::vector<Point>& v);

    virtual ~Polygon() override;

    int verticesCount() const;
    std::vector<Point> getVertices() const;
    bool isConvex() const;

    template<typename... Args>
    Polygon(Args... args);

    template<typename T, typename... Args>
    void build(T value, Args... args);

    template<typename T>
    void build(T value){vertices.push_back(value);};


    virtual void rotate(Point center, double angle) override;
    virtual void reflect(Point center) override;
    virtual void reflect(Line axis) override;
    virtual void scale(Point center, double coefficient) override;

    virtual double perimeter() const override;
    virtual double area() const override;
    virtual bool operator==(const Shape& other) const override;
    virtual bool operator!=(const Shape& other) const override;
    virtual bool isCongruentTo(const Shape& another) const override;
    virtual bool isSimilarTo(const Shape& another) const override;
    virtual bool containsPoint(Point point) const override;

protected:
    std::vector<Point> vertices;
};

class Rectangle : public Polygon {
public:
    Rectangle();
    Rectangle(std::vector<Point> v);
    Rectangle(const Point& a, const Point& b, double c);
    virtual ~Rectangle() override;

    Point center() const;
    std::pair<Line, Line> diagonals() const;
};

class Square : public Rectangle {
public:
    Square(const Point& a, const Point& b);
    virtual ~Square() override;

    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
};

class Triangle : public Polygon {
public:
    Triangle(const Point& a, const Point& b, const Point& c);
    ~Triangle() override;

    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
    Point centroid() const;
    Point orthocenter() const;
    Line EulerLine() const;
    Circle ninePointsCircle() const;

private:
    Line line(uint64_t i) const;
    Line median(uint64_t i) const;
    Line normal(uint64_t i) const;
    Line bisector(uint64_t i) const;
    Point center(uint64_t i) const;
    Point vector(uint64_t i) const;
    Line normal_to_center(uint64_t i) const;
};

void Point::scale(double a) {
    x *= a;
    y *= a;
}

void Point::rotate(double a) {
    a *= pi / 180;
    x = x * cos(a) - y * sin(a);
    y = x * sin(a) + y * cos(a);
}

void Point::reflect(const Point& p) {
    x -= 2 * (x - p.x);
    y -= 2 * (y - p.y);
}

void Point::reflect(const Line& l) { reflect(l.cross(l.normal(*this))); }

Line::Line(const Point& a, const Point& b) {
    if (b.x != a.x) *this = Line(a, (b.y-a.y)/(b.x-a.x));
    else *this = {1,0,-a.x};
}

bool Line::operator==(const Line& another) const {
    bool first = abs((another.a*b) - (another.b*a)) < acc;
    bool second = abs((another.c*b) - (another.b*c)) < acc;
    bool third = abs((another.a*c) - (another.c*a)) < acc;
    return (first && second && third);
}

Line& Line::rotate(const Point& p, const double& alpha) {
    Point pol(0,0);

    if (abs(b) < acc) pol = Point(p.x, (-c -(a * p.x)) / b);
    else Point(p.y, (-c -(b * p.y)) / a);

    pol.x -= p.x;
    pol.y -= p.y;

    pol.rotate(alpha);
    pol.x += p.x; pol.y += p.y;

    *this = Line(p, pol);
    return *this;
}

Point Line::cross(const Line& l) const {
    double d = (a * l.b) - (l.a * b);
    double x = ((b * l.c) - (c * l.b)) / d;
    double y = ((l.a * c) - (l.c * a)) / d;

    return Point(x, y);
}

Line Line::normal(const Point& p) const {
    if (abs(a) > acc) return Line(p, b / a);
    return Line(p, Point(p.x, p.y + 1));
}

Shape::~Shape() {}

std::pair<Line, Line> Ellipse::directrices() const {
    Line x0(f1, f2);
    Point d(f1.x - center().x, f1.y - center().y);

    d.x *= a() / eccentricity() / c();
    d.y *= a() / eccentricity() / c();

    Point left(center().x + d.x, center().y + d.y);
    Point right(center().x - d.x, center().y - d.y);

    return { x0.normal(left), x0.normal(right)};
}

void Ellipse::rotate(Point center, double angle) {
    f1.x -= center.x;
    f1.y -= center.y;

    f2.x -= center.x;
    f2.y -= center.y;

    f1.rotate(angle);
    f2.rotate(angle);

    f1.x += center.x;
    f1.y += center.y;

    f2.x += center.x;
    f2.y += center.y;
}

void Ellipse::reflect(Point center) {
    f1.reflect(center);
    f2.reflect(center);
}

void Ellipse::reflect(Line axis) {
    f1.reflect(axis);
    f2.reflect(axis);
}

void Ellipse::scale(Point center, double coefficient) {
    f1.x = (f1.x - center.x) * coefficient + center.x;
    f1.y = (f1.y - center.y) * coefficient + center.y;

    f2.x = (f2.x - center.x) * coefficient + center.x;
    f2.y = (f2.y - center.y) * coefficient + center.y;

    r *= coefficient;
}

bool Ellipse::operator==(const Shape& other) const {

    const Ellipse* pointer = dynamic_cast<const Ellipse*>(&other);

    if (pointer == nullptr) return false;

    const Ellipse& another = *pointer;

    if (((f1 == another.f1) && (f2 == another.f2) && (abs(r - another.r)<acc)) ||
        ((f1 == another.f2) && (f2 == another.f1) && (abs(r - another.r)<acc))) return true;

    return false;
}

bool Ellipse::operator!=(const Shape& other) const {
    return !(*this == other);
}

bool Ellipse::isCongruentTo(const Shape& other) const {
    const Ellipse* pointer = dynamic_cast<const Ellipse*>(&other);
    if (pointer == nullptr) return false;
    const Ellipse& another = *pointer;
    double dist = sqrt(pow(f1.x-f2.x,2)+pow(f1.y-f2.y,2));
    double anotherDist = sqrt(pow(another.f1.x-another.f2.x,2)+pow(another.f1.y-another.f2.y,2));
    return (abs(dist - anotherDist) < acc) && (abs(r - another.r) < acc);
}

bool Ellipse::isSimilarTo(const Shape& other) const {
    const Ellipse* pointer = dynamic_cast<const Ellipse*>(&other);
    if (pointer == nullptr) return false;
    const Ellipse& another = *pointer;
    double dist = sqrt(pow(f1.x-f2.x,2)+pow(f1.y-f2.y,2));
    double anotherDist = sqrt(pow(another.f1.x-another.f2.x,2)+pow(another.f1.y-another.f2.y,2));
    return (abs(dist/anotherDist) - (r/another.r))<acc;
}

Polygon::Polygon() {}

Polygon::Polygon(const std::vector<Point>& v) : vertices(v) {}

template<typename T, typename... Args>
void Polygon::build(T value, Args... args) {
    vertices.push_back(value);
    build(args...);
}

template<typename... Args>
Polygon::Polygon(Args... args) { build(args...); }

Polygon::~Polygon() {}

int Polygon::verticesCount() const { return vertices.size(); }

std::vector<Point> Polygon::getVertices() const { return vertices; }

bool Polygon::isConvex() const {
    int first = 0; int second = 0;
    for (uint64_t i = 1; i <= vertices.size(); ++i) {
        uint64_t now = (i+vertices.size())%vertices.size();
        uint64_t next = (i+1+vertices.size())%vertices.size();
        uint64_t prev = (i-1+vertices.size())%vertices.size();
        double dx1 = vertices[now].x - vertices[prev].x;
        double dy1 = vertices[now].y - vertices[prev].y;
        double dx2 = vertices[next].x - vertices[now].x;
        double dy2 = vertices[next].y - vertices[now].y;
        if ((dx1*dy2 - dx2*dy1) > acc) ++first;
        else ++second;
    }
    return (first == 0) || (second == 0);
}

void Polygon::rotate(Point center, double angle) {
    for (uint64_t i = 0; i < vertices.size(); ++i) {
        vertices[i].x -= center.x;
        vertices[i].y -= center.y;
        double newvix = vertices[i].x*cos(angle) - vertices[i].y*sin(angle);
        double newviy = vertices[i].x*sin(angle) + vertices[i].y*cos(angle);
        vertices[i].x = newvix + center.x;
        vertices[i].y = newviy + center.y;
    }
}

void Polygon::reflect(Point center) {
    for (Point& p : vertices) { p.reflect(center); }
}

void Polygon::reflect(Line axis) {
    for (Point& p : vertices) { p.reflect(axis); }
}

void Polygon::scale(Point center, double coefficient) {
    for (uint64_t i = 0; i < vertices.size(); ++i) {
        vertices[i].x = (vertices[i].x - center.x)*coefficient + center.x;
        vertices[i].y = (vertices[i].y - center.y)*coefficient + center.y;
    }
}

double Polygon::perimeter() const {
    double res = 0;
    for (uint64_t i = 0; i < vertices.size(); ++i) {
        res += vertices[(i+vertices.size())%vertices.size()].dist(vertices[(i+1+vertices.size())%vertices.size()]);
    }
    return res;
}

double Polygon::area() const {
    double s1 = 0;
    for (uint64_t i = 0; i < (vertices.size()-1); ++i)
        s1 += vertices[i].x*vertices[i+1].y;

    double s2 = 0;
    for (uint64_t i = 0; i < (vertices.size()-1); ++i)
        s2 += vertices[i+1].x*vertices[i].y;

    return 0.5*fabs(s1 + vertices[vertices.size() - 1].x * vertices[0].y - s2 - vertices[0].x * vertices[vertices.size() - 1].y);
}

bool Polygon::operator==(const Shape& other) const {
    const Polygon* pointer = dynamic_cast<const Polygon*>(&other);
    if (pointer == nullptr) return false;
    const Polygon& another = *pointer;
    if (vertices.size() != another.vertices.size()) return false;

    for (uint64_t step = 0; step < vertices.size(); ++step) {
        bool ans = true;
        for (uint64_t index = 0; index < vertices.size(); ++index) {
            if (vertices[index] != another.vertices[(step+index)%vertices.size()]) {
                ans = false;
                break;
            }
        }
        if (ans) return true;
    }
    for (uint64_t step = 0; step < vertices.size(); ++step) {
        bool ans = true;
        for (uint64_t index = 0; index < vertices.size(); ++index) {
            if (vertices[vertices.size()-index-1] != another.vertices[(step+index)%vertices.size()]) {
                ans = false;
                break;
            }
        }
        if (ans) return true;
    }
    return false;
}

bool Polygon::operator!=(const Shape& other) const { return !(*this == other); }

bool Polygon::isCongruentTo(const Shape& other) const {
    return ((area() - other.area()) < acc) && isSimilarTo(other);
}

bool Polygon::isSimilarTo(const Shape& other) const {
    const Polygon* pointer = dynamic_cast<const Polygon*>(&other);

    if (pointer == nullptr) return false;

    const Polygon& another = *pointer;
    const std::vector<Point>& av = another.vertices;

    if (vertices.size() != another.vertices.size()) return false;

    for (uint64_t step = 0; step < vertices.size(); ++step) {
        bool ans = true;
        for (uint64_t index = 0; index < vertices.size(); ++index) {
            int next = (index+1)%vertices.size();
            int nextnext = (index+2)%vertices.size();

            double t1 = (vertices[next].x-vertices[index].x)*(vertices[nextnext].y-vertices[next].y);
            t1 -= (vertices[next].y-vertices[index].y)*(vertices[nextnext].x-vertices[next].x);

            double t2 = (av[(next+step)%av.size()].x-av[(index+step)%av.size()].x)*(av[(nextnext+step)%av.size()].y-av[(next+step)%av.size()].y);
            t2 -= (av[(next+step)%av.size()].y-av[(index+step)%av.size()].y)*(av[(nextnext+step)%av.size()].x-av[(next+step)%av.size()].x);

            if (fabs(t1 - t2) < acc) {
                ans = false;
                break;
            }
        }
        if (ans) return true;
    }
    for (uint64_t step = 0; step < vertices.size(); ++step) {
        bool ans = true;
        for (int64_t index = vertices.size()-1; index > -1; --index) {
            int next = (index-1+vertices.size())%vertices.size();
            int nextnext = (index-2+vertices.size())%vertices.size();

            double t1 = (vertices[next].x-vertices[index].x)*(vertices[nextnext].y-vertices[next].y);
            t1 -= (vertices[next].y-vertices[index].y)*(vertices[nextnext].x-vertices[next].x);

            double t2 = (av[(next+step)%av.size()].x-av[(index+step)%av.size()].x)*(av[(nextnext+step)%av.size()].y-av[(next+step)%av.size()].y);
            t2 -= (av[(next+step)%av.size()].y-av[(index+step)%av.size()].y)*(av[(nextnext+step)%av.size()].x-av[(next+step)%av.size()].x);

            if (fabs(t1 - t2) < acc) {
                ans = false;
                break;
            }
        }
        if (ans) return true;
    }

    return false;
}

bool Polygon::containsPoint(Point point) const {
    double angle = 0;

    for (uint64_t i = 0; i < vertices.size(); ++i) {
        int next = (i + 1) % vertices.size();

        Point pv(vertices[i].x - point.x, vertices[i].y - point.y);
        Point pvnext(vertices[next].x - point.x, vertices[next].y - point.y);

        double v = pv.x * pvnext.y - pv.y * pvnext.x;
        double s = pv.x * pvnext.x + pv.y * pvnext.y;

        angle += atan2(v, s);
    }

    return (fabs(fabs(angle) - 2 * pi) < acc);
}

Rectangle::Rectangle() {}

Rectangle::Rectangle(const Point& a, const Point& b, double c) {
    if (c < 1) c = 1/c;

    double angleAB = (abs(a.x - b.x) < acc) ? (pi / 2) : atan((a.y - b.y)/ ( a.x - b.x));

    double alpha = angleAB + atan(1/c);
    if (alpha > (pi/2)) alpha -= pi;

    Line l1(a, (abs(abs(alpha)-(pi / 2)) < acc) ? Point(a.x, a.y+1) : Point(a.x+(1/tan(alpha)), a.y+tan(alpha)));

    double beta = angleAB + pi - atan(c);
    while (beta > (pi/2)) beta -= pi;

    Line l2(b, (fabs(fabs(beta)-(pi/2))<acc) ? Point(b.x, b.y+1) : Point(b.x+(1/tan(beta)), b.y+tan(beta)));

    Point B = l1.cross(l2);
    Point M = {(a.x+b.x)/2, (a.y+b.y)/2};
    Point D = {2*M.x-B.x, 2*M.y-B.y};

    vertices.push_back(a);
    vertices.push_back(B);
    vertices.push_back(b);
    vertices.push_back(D);
}

Rectangle::~Rectangle() {};

Point Rectangle::center() const {
    return diagonals().first.cross(diagonals().second);
}

std::pair<Line, Line> Rectangle::diagonals() const {
    return {Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3])};
}

Square::Square(const Point& a, const Point& b) {
    double angleAB = (fabs(a.x - b.x) < acc) ? (pi/2) : atan((a.y-b.y)/(a.x-b.x));
    angleAB += pi/2;

    double dist = sqrt(pow(a.x-b.x,2) + pow(a.y-b.y,2))/2;

    Point M = {(a.x+b.x)/2, (a.y+b.y)/2};

    Point B = {cos(angleAB)*dist + M.x, sin(angleAB)*dist + M.y};

    Point D = {2*M.x-B.x, 2*M.y-B.y};

    vertices.push_back(a);
    vertices.push_back(B);
    vertices.push_back(b);
    vertices.push_back(D);
}

Square::~Square() {};

Circle Square::circumscribedCircle() const { return Circle(center(), (vertices[0].dist(vertices[2])) / 2); }

Circle Square::inscribedCircle() const { return Circle(center(), (vertices[0].dist(vertices[1])) /2); }

Triangle::Triangle(const Point& a, const Point& b, const Point& c) : Polygon(a,b,c) {}

Triangle::~Triangle() {}

Line Triangle::normal_to_center(uint64_t i) const {
    uint64_t now = (i+vertices.size())%vertices.size();

    return line(now).normal(center(now));
}

Line Triangle::line(uint64_t i) const {
    uint64_t now = (i+vertices.size())%vertices.size();
    uint64_t next = (i+1+vertices.size())%vertices.size();

    return Line(vertices[now], vertices[next]);
}

Line Triangle::median(uint64_t i) const {
    uint64_t now = (i+vertices.size())%vertices.size();
    uint64_t next = (i+2+vertices.size())%vertices.size();

    return Line(vertices[next], center(now));
}

Line Triangle::normal(uint64_t i) const {
    uint64_t now = (i+vertices.size())%vertices.size();
    uint64_t next = (i+1+vertices.size())%vertices.size();

    return line(next).normal(vertices[now]);
}

Line Triangle::bisector(uint64_t i) const {
    uint64_t now = (i + vertices.size()) % vertices.size();
    uint64_t prev = (i - 1 + vertices.size()) % vertices.size();

    double left_w = vector(now).dist(Point(0,0));
    double right_w = vector(prev).dist(Point(0,0));

    double dx = vertices[now].x + vector(now).x * right_w - left_w*vector(prev).x;
    double dy = vertices[now].y + vector(now).y*right_w - left_w*vector(prev).y;

    return Line(vertices[now], Point(dx, dy));
}

Point Triangle::center(uint64_t i) const {
    uint64_t now = (i + vertices.size()) % vertices.size();
    uint64_t next = (i + 1 + vertices.size()) % vertices.size();

    return Point((vertices[now].x+vertices[next].x)/2, (vertices[now].y+vertices[next].y)/2);
}

Point Triangle::vector(uint64_t i) const {
    uint64_t now = (i + vertices.size()) % vertices.size();
    uint64_t next = (i + 1 + vertices.size()) % vertices.size();

    return Point(vertices[next].x-vertices[now].x, vertices[next].y-vertices[now].y);
}

Circle Triangle::circumscribedCircle() const {
    Point center = normal_to_center(0).cross(normal_to_center(1));
    return Circle(center, center.dist(vertices[0]));
}

Circle Triangle::inscribedCircle() const {
    return Circle(bisector(0).cross(bisector(1)), (2 * area()) / perimeter());
}

Point Triangle::centroid() const { return median(0).cross(median(1)); }

Point Triangle::orthocenter() const { return normal(0).cross(normal(1)); }

Line Triangle::EulerLine() const { return Line(orthocenter(), circumscribedCircle().center());}

Circle Triangle::ninePointsCircle() const {
    Point cross = normal_to_center(0).cross(normal_to_center(1));
    Point center((cross.x + orthocenter().x) / 2 , (cross.y + orthocenter().y) / 2);

    return Circle(center, center.dist(vertices[0]));
}
