#pragma once
#include"util.h"
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
using namespace std;
list<triangle> triangle_cut(vector<point>& points, vector<point>& fa, point near, map<int, vector<point> >& ma_v);


void check_hole(int* is_hole, vector<face>& ans, vector<point>& cur_point, vector<point>& fa, map<pair<int, int>, int>& ma_line, set<int>& is_in_hole);
void dfs(int u, int* head, int* nxt, int* to, int* vis, int* is_hole, vector<int>& hole_s, vector<face>& ans, point* points, vector<point>& fa, map<pair<int, int>, int>& ma_line, set<int>& is_in_hole);
bool check_line(point* points, vector<point>& fa, pair<int, int> tmp_line, pair<int, int> tmp_check_line, map<int, vector<point> >& ma_v);
bool check_line_inter(pair<int, int> cur_line, pair<int, int> pre_line, point* points);
void check_dfs_hole(vector<point>& ma_v, vector<point>& fa, set<face>& s, vector<face>& ans, vector<point>& cur_point, int target);
void check_dfs(int u, point* points, vector<point>& ma_v, vector<point>& fa, set<face>& s, vector<face>& ans, int* head, int* nxt, int* to, int* is_vis, vector<int>& hole_s, int depth, int target);
vector<face> greedy_triangle(vector<float>& ini, vector<point>& fa, vector<int>& point_color);
list<triangle> triangle_cut(vector<point>& points, vector<point>& fa, point near, map<int, vector<point> >& ma_v);
void check_hole(int* is_hole, vector<face>& ans, vector<point>& cur_point, vector<point>& fa, map<pair<int, int>, int>& ma_line, set<int>& is_in_hole);
void dfs(int u, int* head, int* nxt, int* to, int* vis, int* is_hole, vector<int>& hole_s, vector<face>& ans, point* points, vector<point>& fa, map<pair<int, int>, int>& ma_line, set<int>& is_in_hole);
bool check_line(point* points, vector<point>& fa, pair<int, int> tmp_line, pair<int, int> tmp_check_line, map<int, vector<point> >& ma_v);
bool check_line_inter(pair<int, int> cur_line, pair<int, int> pre_line, point* points);
void check_dfs_hole(vector<point>& ma_v, vector<point>& fa, set<face>& s, vector<face>& ans, vector<point>& cur_point, int target);
void check_dfs(int u, point* points, vector<point>& ma_v, vector<point>& fa, set<face>& s, vector<face>& ans, int* head, int* nxt, int* to, int* is_vis, vector<int>& hole_s, int depth, int target);
vector<face> greedy_triangle(vector<float>& ini, vector<point>& fa, vector<int>& point_color);
struct nnangle {
    float angle;
    int index;
    int cur_index;
    bool visible;
    nnangle()
    {
        angle = 0;
        index = -1;
        cur_index = -1;
        visible = false;
    }
};
bool cmp_nnangle(const nnangle& a, const nnangle& b);
struct double_edge {
    point first;
    point second;
    int index;
    double_edge()
    {
        index = -1;
        first = point();
        second = point();
    }
};
vector<point> get_plain_base(vector<point>& points, vector<point>& fa);
bool is_visible(point check_point, point line_first, point line_second, point origin);
vector<face> pcl_greedy_triangle(vector<float>& ini, vector<point>& fa, vector<int>& point_color);


