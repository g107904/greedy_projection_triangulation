#pragma once
#include <iostream>
#include<vector>
#include < fstream >
#include<algorithm>
#include < string >
#include<cstdlib>
#include<stack>
#include<list>
#include<set>
#include<map>
#include<queue>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>
using namespace std;
#define pi 3.1415926
struct point
{
    int pos;
    float indices[3];
    point()
    {
        pos = -1;
        indices[0] = indices[1] = indices[2] = 0;
    }
    point(int pos, float x, float y, float z) { this->pos = pos; this->indices[0] = x; this->indices[1] = y; this->indices[2] = z; }
    point(Eigen::Vector3f in)
    {
        this->pos = -1; this->indices[0] = in[0]; this->indices[1] = in[1]; this->indices[2] = in[2];
    }
    float length()
    {
        float ans = 0;
        for (int i = 0; i < 3; i++)
            ans += this->indices[i] * this->indices[i];
        return sqrt(ans);
    }
    Eigen::Vector3f to_vector()
    {
        Eigen::Vector3f ans;
        for (int i = 0; i < 3; i++)
            ans[i] = this->indices[i];
        return ans;
    }
    float square_length()
    {
        float ans = 0;
        for (int i = 0; i < 3; i++)
            ans += this->indices[i] * this->indices[i];
        return ans;
    }
    point& operator = (const point& a)
    {
        pos = a.pos;
        for (int i = 0; i < 3; i++)
            indices[i] = a.indices[i];
        return *this;
    }
    point operator + (const point& a)
    {
        point ans = point();
        for (int i = 0; i < 3; i++)
            ans.indices[i] = this->indices[i] + a.indices[i];
        return ans;
    }
    point operator / (int num)
    {
        point ans = point();
        for (int i = 0; i < 3; i++)
            ans.indices[i] = this->indices[i] / (float)num;
        return ans;
    }

    point operator - (const point& a)
    {
        point ans = point();
        for (int i = 0; i < 3; i++)
            ans.indices[i] = this->indices[i] - a.indices[i];
        return ans;
    }
    point operator - ()
    {
        return point() - *this;
    }
    point operator * (float a)
    {
        point ans = point();
        for (int i = 0; i < 3; i++)
            ans.indices[i] = this->indices[i] * a;
        return ans;
    }
    bool operator == (const point& a)
    {
        bool flag = true;
        if (this->pos != a.pos)
            flag = false;
        for (int i = 0; i < 3 && flag; i++)
            if (this->indices[i] != a.indices[i])
                flag = false;
        return flag;
    }
    point norm()
    {
        float sum = 0;
        for (int i = 0; i < 3; i++)
            sum += indices[i] * indices[i];
        for (int i = 0; i < 3; i++)
            indices[i] = indices[i] / sqrt(sum);
        return *this;
    }
};
bool cmpx(const point& a, const point& b);
bool cmpy(const point& a, const point& b);
bool cmpz(const point& a, const point& b);
bool cmppos(const point& a, const point& b);
bool cmpxsuby(const point& a, const point& b);
bool cmpxandy(const point& a, const point& b);
//typedef bool (*cmp)(const point& a, const point& b);
//cmp ccmp[3] = { cmpx,cmpy,cmpz };
struct line
{
    point* l;
    point* r;
    line() { l = nullptr; r = nullptr; }
    line(point* l, point* r)
    {
        this->l = l;
        this->r = r;
    }
};
struct node
{
    int depth;
    node* left;
    node* right;
    point cur;
    node(int depth, point cur)
    {
        this->depth = depth;
        this->cur = cur;
        this->left = nullptr;
        this->right = nullptr;
    }
};
struct triangle
{
    point points[3];
    triangle()
    {

    }
    triangle(const point& a, const point& b, const point& c)
    {
        this->points[0] = a; this->points[1] = b; this->points[2] = c;
    }

    inline bool operator == (const triangle& ps) const
    {
        bool flag = true;
        set<int> s;
        for (int i = 0; i < 3; i++)
            s.insert(ps.points[i].pos);
        for (int i = 0; i < 3 && flag; i++)
            flag = s.find(this->points[i].pos) != s.end();
        s.clear();
        return flag;
    }
};
struct face
{
    int indices[3];
    face()
    {
        for (int i = 0; i < 3; i++)
            indices[i] = 0;
    }
    face(int x, int y, int z)
    {
        indices[0] = x; indices[1] = y; indices[2] = z;
        sort(indices, indices + 3, less<int>());
    }
    face(triangle* x)
    {
        for (int i = 0; i < 3; i++)
            indices[i] = x->points[i].pos;
        sort(indices, indices + 3, less<int>());
    }

    friend bool operator < (const face& a, const face& b)
    {
        for (int i = 0; i < 3; i++)
            if (a.indices[i] != b.indices[i])
                return a.indices[i] < b.indices[i];
        return false;
    }
};
bool cmpind(const face& a, const face& b);
struct pqcmp {
    template<typename T, typename U>
    bool operator()(T const& left, U const& right) {
        if (left.first < right.first) return true;
        return false;
    }
};
struct pq_node {
    float dis;
    int pos;
    pq_node()
    {
        dis = 0; pos = -1;
    }
    pq_node(float dis, int pos)
    {
        this->dis = dis;
        this->pos = pos;
    }
    friend bool operator > (const pq_node& a, const pq_node& b)
    {
        return a.dis > b.dis;
    }
    friend bool operator < (const pq_node& a, const pq_node& b)
    {
        return a.dis < b.dis;
    }
};
node* kd_build(point* points, int len, int depth);
void kd_add(node* curnode, point target);
float l2_distance(const point& a, const point& b);
float x_mul(point a, point b);
point x_mul3d(point t1, point t2);
point* kd_nearst_search(node* curnode, point target);
vector<point> kd_k_search(node* curnode, point target, int maxnum);
bool is_right_point(point a, point b, point c, vector<point>& fa);
float p_distance(point a, line* b);
float get_angle(line* a, line* b);
float get_angle(point* a, point* b, point* c);
float get_min_angle(point* a, point* b, point* c);
bool blur_triangle(triangle* t, float mmin, float mmax);
bool is_in_triangle(triangle* tr, point* pt);
bool is_on_triangle(triangle* tr, point* pt);

bool cmp_ld(const point& a, const point& b, point ld_point);
point get_fa(vector<point>& points);
float point_to_plain(point* cur, point* fa, float bias);
bool is_in_circle(point p, triangle t);
void get_base(point fa, point* x_base, point* y_base);
point get_plain_point(point plain_fa, float d, point p, point fa);
float dot_mul(point a, point b);
bool cmp_zero(point a, point b, point zero);
int check_if_triangle(point* points, int pre, int i, int late, int terminal, int len);
int find_pa(vector<int>& pa, int pos);
void smallest_tree(vector<point>& fa, point* points, node* root);
float get_angle_2d(line* a, line* b);
