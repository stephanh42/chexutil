#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>

#include <string>
#include <vector>
#include <cctype>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <memory>
#include <unordered_set>
#include <map>

namespace py = pybind11;

struct InvalidHex : public std::invalid_argument
{
  InvalidHex(const char *what) : std::invalid_argument(what) {}
};

struct DivMod
{
  int div;
  int mod;

  DivMod(int p, int q)
    : div(p/q), mod(p%q)
  {
    if (mod < 0) {
      mod += q;
      div -= 1;
    }
  }
};

std::pair<int,int> tiled_range(int lo, int hi, int tile_size)
{
  return std::pair<int,int>(DivMod(lo, tile_size).div, DivMod(hi + tile_size - 1, tile_size).div);
}

std::pair<int,int> make_range(int x, int width, int bloat, int grid_size)
{
  return tiled_range(x + grid_size - 1 - bloat, x + width + bloat, grid_size);
}

typedef std::function<int(int)> RandomFunction;

RandomFunction make_random_function(py::object random)
{
  if (random.is_none()) {
    return RandomFunction([](int n) { return rand() % n; });
  } else {
    py::function randrange = random.attr("randrange").cast<py::function>();
    return RandomFunction([randrange](int n) { return randrange(0, py::cast(n)).cast<int>(); });
  }
}

class Hex {
  public:
    Hex() : x(0), y(0) {}
    Hex(int the_x, int the_y) : x(the_x), y(the_y) {}

    int x;
    int y;
    static std::array<Hex, 6> orig_neighbours;

    static Hex create_checked(int x, int y)
    {
      if ((x ^ y) & 1) {
        throw InvalidHex("Sum of x and y must be even.");
      }
      return Hex(x, y);
    }

    bool operator==(const Hex &other) const
    { 
      return (x == other.x) && (y == other.y);
    }

    Hex operator+(const Hex &other) const
    {
      return Hex(x + other.x, y + other.y);
    }

    Hex operator-(const Hex &other) const
    {
      return Hex(x - other.x, y - other.y);
    }

    Hex operator-() const
    {
      return Hex(-x, -y);
    }

    py::list neighbours() const
    {
      py::list result;
      for (const Hex &nb: orig_neighbours) {
        result.append((*this) + nb);
      }
      return result;
    }

    // Return a random neighbour of this hexagon.
    Hex random_neighbour(const RandomFunction &random) const
    {
      return (*this) + orig_neighbours[random(6)];
    }

    // Returns a list of length N+1 since it includes the start point.
    py::list random_walk(int N, const RandomFunction &random) const
    {
      py::list result;
      Hex position = *this;
      result.append(py::cast(position));
      for (int i = 0; i < N; i++) {
        position = position.random_neighbour(random);
        result.append(py::cast(position));
      }
      return result;
    }

    std::string repr() const
    {
      std::ostringstream out;
      out << "Hex(" << x << ", " << y << ")";
      return out.str();
    }

    // Distance in number of hexagon steps.
    // Direct neighbours of this hex have distance 1.
    int distance(const Hex &other) const
    {
      int dx = std::abs(x - other.x);
      int dy = std::abs(y - other.y);
      return dy + std::max(0, (dx - dy)>>1);
    }

    //Given a hex return the hex when rotated 60° counter-clock-wise around the origin.
    Hex rotate_left() const
    {
        return Hex((x - 3 * y) >> 1, (x + y) >> 1);
    }

    // Given a hex return the hex when rotated 60° clock-wise around the origin.
    Hex rotate_right() const
    {
      return Hex((x + 3 * y) >> 1, (y - x) >> 1);
    }

};

std::array<Hex, 6> Hex::orig_neighbours{{Hex{2, 0}, Hex{1, 1}, Hex{-1, 1}, Hex{-2, 0}, Hex{-1, -1}, Hex{1, -1}}};

typedef Hex (*HexTransform)(const Hex&);

const std::array<HexTransform, 6> rotations { {
  [](const Hex &hex) { return hex; },
  [](const Hex &hex) { return hex.rotate_left(); },
  [](const Hex &hex) { return -hex.rotate_right(); },
  [](const Hex &hex) { return -hex; },
  [](const Hex &hex) { return -hex.rotate_left(); },
  [](const Hex &hex) { return hex.rotate_right(); },
} };

struct Rectangle {
  Rectangle(int ax, int ay, int awidth, int aheight)
  : x(ax), y(ay), width(awidth), height(aheight)
  {
    if (empty()) {
      x = y = width = height = 0;
    }
  }

  int x;
  int y;
  int width;
  int height;

  static Rectangle from_coordinates(int x1, int y1, int x2, int y2)
  {
    return Rectangle(x1, y1, x2-x1, y2-y1);
  }

  bool operator==(const Rectangle &other) const
  {
    return (x == other.x) && (y == other.y) && (width == other.width) && (height == other.height);
  }

  int x2() const
  {
    return x + width;
  }

  int y2() const
  {
    return y + height;
  }

  bool empty() const
  {
    return (width <= 0) || (height <= 0);
  }

  Rectangle translated(int dx, int dy) const
  {
    return Rectangle(x+dx, y+dy, width, height);
  }

  Rectangle operator|(const Rectangle &other) const
  {
    if (empty()) {
      return other;
    } else if (other.empty()) {
      return *this;
    } else {
      return from_coordinates(
        std::min(x, other.x),
        std::min(y, other.y),
        std::max(x2(), other.x2()),
        std::max(y2(), other.y2()));
    }
  }

  Rectangle operator&(const Rectangle &other) const
  {
    return from_coordinates(
        std::max(x, other.x),
        std::max(y, other.y),
        std::min(x2(), other.x2()),
        std::min(y2(), other.y2()));
  }

  std::string repr() const
  {
    std::ostringstream out;
    out << "Rectangle(" << x << ", " << y << ", " << width << ", " << height << ")";
    return out.str();
  }
};

struct HexIterator
{
  typedef int Sentinel;

  int x, y, x1, x2;

  Hex operator*() { return Hex(x, y); }

  bool operator==(Sentinel y2) const
  {
    return y == y2;
  }

  HexIterator &operator++()
  {
    x += 2;
    if (x >= x2) {
      x = x1;
      y++;
      if ((x + y) & 1) {
        if ((x + 1) == x2) {
          y++;
        } else {
          x++;
        }
      }
    }
    return *this;
  }
};

struct HexRange
{
  HexRange(int ax1, int ax2, int ay1, int ay2)
  : x1(ax1), x2(ax2), y1(ay1), y2(ay2)
  {
    if ((x1 + 1) == x2) {
      y1 += (x1 + y1) & 1;
      y2 += (x1 + y2) & 1;
    } else if ((y1 + 1) == y2) {
      x1 += (x1 + y1) & 1;
      x2 += (x2 + y1) & 1;
    }
    if ((x1 >= x2) || (y1 >= y2)) {
      x1 = x2 = y1 = y2 = 0;
    }
  }

  int x1, x2, y1, y2;

  py::iterator iter() const
  {
    int x = x1 + ((x1 + y1) & 1);
    return py::make_iterator<py::return_value_policy::copy, HexIterator, HexIterator::Sentinel>(
        HexIterator{x, y1, x1, x2}, y2);
  }

  bool contains(const Hex &hex) const
  {
    return (x1 <= hex.x) && (hex.x < x2) && (y1 <= hex.y) && (hex.y < y2);
  }

  std::string repr() const
  {
    std::ostringstream out;
    out << "HexRange(x1=" << x1 << ", x2=" << x2 << ", y1=" << y1 << ", y2=" << y2 << ")";
    return out.str();
  }
};

std::array<Hex, 6> origin_corners = {{Hex(1, 1), Hex(0, 2), Hex(-1, 1), Hex(-1, -1), Hex(0, -2), Hex(1, -1)}};

class HexGrid {
    static const double hex_factor;

  public:
    HexGrid(int w, int h) : width(w), height(h) {}
    HexGrid(int w) : width(w), height(static_cast<int>(std::round(w * hex_factor))) {}

    int width;
    int height;

    // Get the center (as (x, y) tuple) of a hexagon.
    std::pair<int, int> center(const Hex &hex) const
    {
        return std::pair<int,int>(hex.x*width, 3*height*hex.y);
    }

    // Get the bounding box (as a Rectangle) of a hexagon.
    Rectangle bounding_box(const Hex &hex) const
    {
        std::pair<int, int> c = center(hex);
        return Rectangle(c.first - width, c.second - 2*height, 2*width, 4*height);
    }

    py::list corners(const Hex &hex) const
    {
      py::list result;
      int x0 = hex.x;
      int y0 = 3 * hex.y;
      for (const Hex &corner: origin_corners) {
	result.append(py::cast(std::pair<int, int>(width * (corner.x + x0), height * (corner.y + y0))));
      }
      return result;
    }

    // Given pixel coordinates x and y, get the hexagon under it.
    Hex hex_at_coordinate(int x, int y)
    {
      DivMod x0(x, width);
      DivMod y0(y, 3 * height);

      if (((x0.div + y0.div) & 1) == 0) {
        if (width * y0.mod < height * (2 * width - x0.mod)) {
          return Hex(x0.div, y0.div);
        } else {
          return Hex(x0.div + 1, y0.div + 1);
        }
      } else if (width * y0.mod < height * (width + x0.mod)) {
        return Hex(x0.div + 1, y0.div);
      } else {
        return Hex(x0.div, y0.div + 1);
      }
    }

    // Return a sequence with the hex coordinates in the rectangle."""
    HexRange hexes_in_rectangle(const Rectangle &rectangle)
    {
      if (rectangle.empty()) {
        return HexRange{0, 0, 0, 0};
      } else {
        std::pair<int,int> x_range = make_range(rectangle.x, rectangle.width, width, width);
        std::pair<int,int> y_range = make_range(rectangle.y, rectangle.height, 2*height, 3*height);
        return HexRange{x_range.first, x_range.second, y_range.first, y_range.second};
      }
    }
 
    std::string repr() const
    {
      std::ostringstream out;
      out << "HexGrid(" << width << ", " << height << ")";
      return out.str();
    }
};

const double HexGrid::hex_factor = std::sqrt(1.0/3.0);

py::tuple make_rotations()
{
  py::list result;
  for (auto rotation : rotations) {
    result.append(py::cpp_function(rotation, py::arg("hex")));
  }
  return py::tuple(result);
}

class FovTree 
{
  private:
    static const std::array<Hex, 4> corners;
    static const std::array<Hex, 3> neighbours;

    bool m_has_cached_successors;
    std::vector<FovTree> m_cached_successors;

    Hex m_hexagon;
    double m_angle1;
    double m_angle2;
    unsigned m_direction;

    std::array<Hex, 6> m_hexagons;
    int m_distance;

  public:
    FovTree(const Hex &hexagon, unsigned direction, double angle1, double angle2)
    : m_has_cached_successors(false),
      m_hexagon(hexagon), m_angle1(angle1), m_angle2(angle2), m_direction(direction),
      m_distance(hexagon.distance(Hex()))
    {
      size_t pos = 0;
      for (auto rotation : rotations) {
        m_hexagons[pos] = rotation(hexagon);
        pos++;
      }
    }

    double get_angle(const Hex &corner) const
    {
      return (3*m_hexagon.y + corner.y)/double(m_hexagon.x + corner.x);
    }

    static const unsigned all_directions;

    void field_of_view(const Hex &offset, unsigned direction, const py::function &transparent, int max_distance, 
        py::dict &visible, py::function &visible_get)
    {
      if (m_distance > max_distance) {
        return;
      }
      py::object hexagon = py::cast(offset + m_hexagons[direction]);
      if (PyObject_IsTrue(transparent(hexagon).ptr())) {
        visible[hexagon] = all_directions;
        for (FovTree &succ : successors()) {
          succ.field_of_view(offset, direction, transparent, max_distance, visible, visible_get);
        }
      } else {
        int directions = 1 << ((m_direction + direction) % 6);
        visible[hexagon] = directions | visible_get(hexagon, 0).cast<int>();
      }
    }

    std::vector<FovTree> &successors()
    {
      if (!m_has_cached_successors) {
        m_has_cached_successors = true;
        std::array<double, 4> angles;
        for (unsigned i = 0; i < 4; i++) {
          angles[i] = get_angle(corners[i]);
        }
        for (unsigned i = 0; i < 3; i++) {
          double c1 = std::max(m_angle1, angles[i]);
          double c2 = std::min(m_angle2, angles[i+1]);
          if (c1 < c2) {
            const Hex &nb = neighbours[i];
            m_cached_successors.push_back(FovTree(m_hexagon + nb, (i+5) % 6, c1, c2));
          }
        }
      }
      return m_cached_successors;
    }
};

const unsigned FovTree::all_directions = (1 << 6) - 1;
const std::array<Hex, 4> FovTree::corners = {{Hex(0, -2), Hex(1, -1), Hex(1, 1), Hex(0, 2)}};
const std::array<Hex, 3> FovTree::neighbours = {{Hex(1, -1), Hex(2, 0), Hex(1, 1)}};

const char *field_of_view_doc = R"(
Calculate field-of-view.

transparent  -- from a Hex to a boolean, indicating of the Hex is transparent
max_distance -- maximum distance you can view
visible      -- if provided, should be a dict which will be filled and returned

Returns a dict which has as its keys the hexagons which are visible.
The value is a bitmask which indicates which sides of the hexagon are visible.
The bitmask is useful if you want to use this function also to compute light sources.

view_set = player_pos.field_of_view(...)
light_set = light_source.field_of_view(...)

# Is pos visible?
if view_set.get(pos, 0) & light_set.get(pos, 0):
    # yes it is
)";

py::dict field_of_view(const Hex &hex, py::function transparent, int max_distance, py::object visible)
{
  static FovTree fovtree(Hex(2, 0), 0, -1.0, 1.0);

  py::dict visible_dict = visible.is_none() ? py::dict() : visible.cast<py::dict>();
  visible_dict[py::cast(hex)] = FovTree::all_directions;
  py::function visible_get = visible_dict.attr("get").cast<py::function>();
  for (unsigned direction = 0; direction < 6; direction++) {
    fovtree.field_of_view(hex, direction, transparent, max_distance, visible_dict, visible_get);
  }
  return visible_dict;
}

// custom specialization of std::hash can be injected in namespace std
namespace std
{
  template<> struct hash<Hex>
  {
    typedef Hex argument_type;
    typedef std::size_t result_type;
    result_type operator()(const argument_type& hex) const
    {
      return (hex.x >> 1) + (hex.y * 2738050981);
    }
  };

  template<> struct hash<Rectangle>
  {
    typedef Rectangle argument_type;
    typedef std::size_t result_type;
    result_type operator()(const argument_type& rect) const
    {
      return rect.x + rect.y * 414353063 + rect.width * 371838251 + rect.height * 750515191;
    }
  };
}

const char *HexPathFinder_doc =
R"(A* path-finding on the hex grid.
All positions are represented as Hex objects.

Important data attributes: 
found -- True if path-finding is complete and we found a path
done  -- True if path-finding is complete: we either found a path or know there isn't one
path  -- The path, as a tuple of positions from start to destination (including both). Empty tuple if found is False.
)";

struct HexPath
{
  typedef std::shared_ptr<HexPath> Ptr;

  HexPath(const Hex &pos, Ptr the_next)
    : position(pos), next(the_next)
  {}

  Hex position;
  Ptr next;
};

class HexPathFinder {
  public:
    bool found;
    bool done;
    std::vector<Hex> path;

    typedef std::function<bool(Hex)> PassableFunction;
    typedef std::function<double(Hex)> CostFunction;

  private:
    Hex m_start;
    Hex m_destination;
    py::function m_passable;
    CostFunction m_cost;

    struct OpenItem {
      double cost_so_far;
      Hex position;
      HexPath::Ptr path;
    };

    std::unordered_set<Hex> m_closed;
    std::map<double, std::vector<OpenItem> > m_open;

    void add_to_open_set(const OpenItem &item)
    {
      double h = heuristic(item.position) + item.cost_so_far;
      m_open[h].push_back(item);
    }

    bool pop_from_open_set(OpenItem &item)
    {
      while (true) {
        auto begin_iter = m_open.begin();
        if (begin_iter == m_open.end()) {
          return false;
        } else {
          std::vector<OpenItem> &vec = begin_iter->second;
          if (vec.empty()) {
            m_open.erase(begin_iter);
            continue;
          } else {
            item = vec.back();
            vec.pop_back();
            return true;
          }
        }
      }
    }

    double heuristic(const Hex &position)
    {
      return m_destination.distance(position);
    }

  public:
    static CostFunction makeCostFunction(py::object cost)
    {
      if (cost.is_none()) {
        return [](Hex) { return 1.0; };
      } else {
        return [cost](Hex hex) { return cost(py::cast(hex)).cast<double>(); };
      }
    }

    bool is_passable(const Hex &hex)
    {
      return PyObject_IsTrue(m_passable(py::cast(hex)).ptr());
    }

    HexPathFinder(const Hex &start, const Hex &destination, py::function passable, py::object cost)
    : found(false), done(false),
      m_start(start),
      m_destination(destination),
      m_passable(passable),
      m_cost(makeCostFunction(cost))
    {
      add_to_open_set(OpenItem{0.0, start, nullptr});
    }

    std::vector<Hex> compute_path(HexPath::Ptr path) const
    {
      std::vector<Hex> result;
      while (path) {
        result.push_back(path->position);
        path = path->next;
      }
      std::reverse(result.begin(), result.end());
      return result;
    }

    void run_n(int n)
    {
      // Run at most n path-finding steps
      // This method does a bounded amount of work, and is therefore useful
      // if pathfinding must be interleaved with interactive behaviour or may
      // be interrupted.
      for (int i = 0; i < n; i++) {
        OpenItem current;
        if (!pop_from_open_set(current)) {
          done = true;
          return;
        }
        if (m_closed.find(current.position) != m_closed.end()) {
          continue;
        }
        HexPath::Ptr new_path = std::make_shared<HexPath>(HexPath(current.position, current.path));

        if (current.position == m_destination) {
          path = compute_path(new_path);
          found = done = true;
          m_open.clear();
          return;
        }
        m_closed.insert(current.position);

        for (const Hex &delta : Hex::orig_neighbours) {
          Hex new_pos = current.position + delta;
          if (!is_passable(new_pos) || (m_closed.find(new_pos) != m_closed.end())) {
            continue;
          }
          double new_cost = current.cost_so_far + m_cost(new_pos);
          add_to_open_set(OpenItem{new_cost, new_pos, new_path});
        }
      }
    }

    void run()
    {
      // Run path-finding until done, that is, we either found a path or know there isn't one.
      while (!done) {
        run_n(1000);
      }
    }
};

const char *find_path_doc =
R"(Perform path-finding.
self        -- Starting position for path finding.
destination -- Destination position for path finding.
passable    -- Function of one position, returning True if we can move through this hex.
cost        -- cost function for moving through a hex. Should return a value ≥ 1. By default all costs are 1.
)";

std::vector<Hex> 
find_path(const Hex &start, const Hex &destination, py::function passable, py::object cost)
{
  HexPathFinder pathFinder(start, destination, passable, cost);
  pathFinder.run();
  return pathFinder.path;
}

PYBIND11_MODULE(chexutil, m) {
  m.doc() = R"(
Classes and functions to deal with hexagonal grids.

This module assumes that the hexagonal grid is aligned with the x-axis.
If you need it to be aligned with the y-axis instead, you will have to
swap x and y coordinates everywhere.
)";

  py::register_exception<InvalidHex>(m, "InvalidHex");
  py::class_<Hex> hexClass =
    py::class_<Hex>(m, "Hex",
    "A single hexagon in a hexagonal grid.")
    .def(py::init(&Hex::create_checked), py::arg("x"), py::arg("y"))
    .def(py::init<>())
    .def_readonly("x", &Hex::x)
    .def_readonly("y", &Hex::y)
    .def(py::self == py::self)
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(- py::self)
    .def(hash(py::self))
    .def("__repr__", &Hex::repr)
    .def("neighbours", &Hex::neighbours)
    .def("distance", &Hex::distance, py::arg("other_hex"),
        R"(Distance in number of hexagon steps.
      Direct neighbours of this hex have distance 1.)")
    .def("rotate_left", &Hex::rotate_left,
        "Given a hex return the hex when rotated 60° counter-clock-wise around the origin.")
    .def("rotate_right", &Hex::rotate_right,
        "Given a hex return the hex when rotated 60° clock-wise around the origin.")
    .def("random_neighbour", 
        [](const Hex &hex, py::object random) { return hex.random_neighbour(make_random_function(random)); },
        py::arg("random")=py::none(),
        "Return a random neighbour of this hexagon.")
    .def("random_walk", 
        [](const Hex &hex, int N, py::object random) { return hex.random_walk(N, make_random_function(random)); },
        py::arg("N"), py::arg("random")=py::none(),
        "Returns a list of length N+1 since it includes the start point.")
    .def("__iter__", [](const Hex &hex) {
      std::array<int, 2> coords{{hex.x, hex.y}};
      return py::iter(py::cast(coords));
    })
    .def("field_of_view", &field_of_view, py::arg("transparent"), py::arg("max_distance"), py::arg("visible")=py::none(),
        field_of_view_doc)
    .def("find_path", &find_path, py::arg("destination"), py::arg("passable"), py::arg("cost")=py::none(),
        find_path_doc)
    .def(py::pickle(
	  [](const Hex &hex) { // __getstate__
  	    /* Return a tuple that fully encodes the state of the object */
	    return py::make_tuple(hex.x, hex.y);
	  },
	  [](py::tuple t) { // __setstate__
	    if (t.size() != 2) throw std::runtime_error("Invalid state!");
	    return Hex(t[0].cast<int>(), t[1].cast<int>());
	  }))
    ;
  hexClass.attr("rotations") = make_rotations();

  m.attr("origin") = py::cast(Hex());

  py::class_<Rectangle>(m, "Rectangle",
    R"(Represents a rectangle.
    x, y   -- position of lower-left corner
    width  -- width of rectangle
    height -- height of rectangle
    )")
    .def(py::init<int, int, int, int>(), py::arg("x"), py::arg("y"), py::arg("width"), py::arg("height"))
    .def(py::init([](int w, int h) { return Rectangle(0, 0, w, h); }), py::arg("width"), py::arg("height"))
    .def_static("from_coordinates", &Rectangle::from_coordinates)
    .def_readonly("x", &Rectangle::x)
    .def_readonly("y", &Rectangle::y)
    .def_readonly("width", &Rectangle::width)
    .def_readonly("height", &Rectangle::height)
    .def_property_readonly("x2", &Rectangle::x2)
    .def_property_readonly("y2", &Rectangle::y2)
    .def(py::self == py::self)
    .def(hash(py::self))
    .def("empty", &Rectangle::empty, "Check if Rectangle is empty.")
    .def("translated", &Rectangle::translated, py::arg("dx"), py::arg("dy"))
    .def(py::self | py::self)
    .def(py::self & py::self)
    .def("__repr__", &Rectangle::repr)
    .def("__iter__", [](const Rectangle &rect) {
      std::array<int, 4> coords{{rect.x, rect.y, rect.width, rect.height}};
      return py::iter(py::cast(coords));
    })
    .def(py::pickle(
	  [](const Rectangle &rect) { // __getstate__
  	    /* Return a tuple that fully encodes the state of the object */
	    return py::make_tuple(rect.x, rect.y, rect.width, rect.height);
	  },
	  [](py::tuple t) { // __setstate__
	    if (t.size() != 4) throw std::runtime_error("Invalid state!");
	    return Rectangle(t[0].cast<int>(), t[1].cast<int>(), t[2].cast<int>(), t[3].cast<int>());
	  }))
    ;

  py::class_<HexPathFinder>(m, "HexPathFinder", HexPathFinder_doc)
    .def(py::init<const Hex&, const Hex&, py::function, py::object>(),
        py::arg("start"), py::arg("destination"), py::arg("passable"), py::arg("cost")=py::none(),
        R"(Create a new HexPathFinder object.
        start       -- Starting position for path finding.
        destination -- Destination position for path finding.
        passable    -- Function of one position, returning True if we can move through this hex.
        cost        -- cost function for moving through a hex. Should return a value ≥ 1. By default all costs are 1.
        )")
    .def_readonly("done", &HexPathFinder::done)
    .def_readonly("found", &HexPathFinder::found)
    .def_readonly("path", &HexPathFinder::path)
    .def("run", &HexPathFinder::run,
        "Run path-finding until done, that is, we either found a path or know there isn't one.")
    .def("run_n", &HexPathFinder::run_n, py::arg("n"),
        R"(Run at most n path-finding steps.
        This method does a bounded amount of work, and is therefore useful
        if pathfinding must be interleaved with interactive behaviour or may
        be interrupted.
        )")
     ;

  py::class_<HexRange>(m, "HexRange",
      R"(A 2D range (rectangle) of hexes.)")
    .def(py::init([](int x1, int x2, int y1, int y2) { return HexRange{x1, x2, y1, y2}; }), 
         py::arg("x1"), py::arg("x2"), py::arg("y1"), py::arg("y2"))
    .def_readonly("x1", &HexRange::x1)
    .def_readonly("x2", &HexRange::x2)
    .def_readonly("y1", &HexRange::y1)
    .def_readonly("y2", &HexRange::y2)
    .def("__iter__", &HexRange::iter)
    .def("__repr__", &HexRange::repr)
    .def("__contains__", &HexRange::contains)
    ;

  py::class_<HexGrid>(m, "HexGrid",
    R"(Represents the dimensions of a hex grid as painted on the screen.
    The hex grid is assumed to be aligned horizontally, like so:
       / \ / \ / \ 
      |   |   |   |
       \ / \ / \ /
    The center of hex (0, 0) is assumed to be on pixel (0, 0).
    The hexgrid is determined by width and height, which are the screen coordinates
    of the upper-right corner of the central hex.

    To have equilateral hexes, width:height should be approximately √3 : 1.
    If you only pass in width to the constructor, the height is computed to be
    an integer as close as possible to width / √3 .
    )")
    .def(py::init<int>())
    .def(py::init<int, int>())
    .def("__repr__", &HexGrid::repr)
    .def_readonly("width", &HexGrid::width)
    .def_readonly("height", &HexGrid::height)
    .def("center", &HexGrid::center,
        "Get the center (as (x, y) tuple) of a hexagon.")
    .def("corners", &HexGrid::corners,
        "Get the 6 corners (in pixel coordinates) of the hex.")
    .def("bounding_box", &HexGrid::bounding_box,
        "Get the bounding box (as a Rectangle) of a hexagon.")
    .def("hex_at_coordinate", &HexGrid::hex_at_coordinate, py::arg("x"), py::arg("y"),
        "Given pixel coordinates x and y, get the hexagon under it.")
    .def("hexes_in_rectangle", &HexGrid::hexes_in_rectangle, py::arg("rectangle"),
        "Return a sequence with the hex coordinates in the rectangle.")
    ;
 }
