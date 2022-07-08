#include "util.h"
bool cmpx(const point& a, const point& b) { return a.indices[0] < b.indices[0]; }
bool cmpy(const point& a, const point& b) { return a.indices[1] < b.indices[1]; }
bool cmpz(const point& a, const point& b) { return a.indices[2] < b.indices[2]; }
bool cmppos(const point& a, const point& b) { return a.pos < b.pos; }
bool cmpxsuby(const point& a, const point& b) { return a.indices[0] - a.indices[1] < b.indices[0] - b.indices[1]; }
bool cmpxandy(const point& a, const point& b) { return a.indices[0] + a.indices[1] < b.indices[0] + b.indices[1]; }
typedef bool (*cmp)(const point& a, const point& b);
cmp ccmp[3] = { cmpx,cmpy,cmpz };
bool cmpind(const face& a, const face& b)
{
    if (a.indices[0] != b.indices[0])
        return a.indices[0] < b.indices[0];
    if (a.indices[1] != b.indices[1])
        return a.indices[1] < b.indices[1];
    return a.indices[2] < b.indices[2];
}
node* kd_build(point* points, int len, int depth)
{

    if (points == nullptr || len <= 0)
        return nullptr;
    int cut_dim = depth % 3;
    int mid = len / 2;
    sort(points, points + len, ccmp[cut_dim]);
    node* cur = new node(depth, *(points + mid));
    cur->left = kd_build(points, mid, depth + 1);
    cur->right = kd_build(points + mid + 1, len - mid - 1, depth + 1);
    return cur;

}
void kd_add(node* curnode, point target)
{
    stack<node*> s;
    while (curnode != nullptr)
    {
        s.push(curnode);
        if (target.indices[curnode->depth % 3] <= (curnode->cur).indices[curnode->depth % 3])
            curnode = curnode->left;
        else
            curnode = curnode->right;
    }
    curnode = s.top();
    if (target.indices[curnode->depth % 3] <= (curnode->cur).indices[curnode->depth % 3])
        curnode->left = kd_build(&target, 1, curnode->depth + 1);
    else
        curnode->right = kd_build(&target, 1, curnode->depth + 1);
}
float l2_distance(const point& a, const point& b)
{
    float s = 0;
    for (int i = 0; i < 3; i++)
        s += (a.indices[i] - b.indices[i]) * (a.indices[i] - b.indices[i]);
    return sqrt(s);
}
float x_mul(point a, point b)
{
    return a.indices[0] * b.indices[1] - a.indices[1] * b.indices[0];
}
point x_mul3d(point t1, point t2)
{
    point ans;
    ans = point(-1, t1.indices[1] * t2.indices[2] - t1.indices[2] * t2.indices[1], t1.indices[2] * t2.indices[0] - t1.indices[0] * t2.indices[2], t1.indices[0] * t2.indices[1] - t1.indices[1] * t2.indices[0]);
    return ans;
}
point* kd_nearst_search(node* curnode, point target)
{

    stack<node*> s;
    while (curnode != nullptr)
    {
        s.push(curnode);
        if (target.indices[curnode->depth % 3] <= (curnode->cur).indices[curnode->depth % 3])
            curnode = curnode->left;
        else
            curnode = curnode->right;

    }
    point* nearst = new point();
    /*
    nearst = &(s.top())->cur;
    s.pop();
    float dist = l2_distance(*nearst, target);
    if (nearst->pos == target.pos)
    */

    float dist = 1000.0;
    set<int> is_in;
    while (!s.empty())
    {
        node* tmp = s.top(); s.pop();
        if (tmp->left == nullptr && tmp->right == nullptr)
        {
            if (dist > l2_distance(tmp->cur, target) && tmp->cur.pos != target.pos)
            {
                nearst = &(tmp->cur);
                dist = l2_distance(tmp->cur, target);

            }
        }
        else
        {
            int cut_dim = tmp->depth % 3;
            if (fabs(tmp->cur.indices[cut_dim] - target.indices[cut_dim]) < dist)
            {
                if (dist > l2_distance(tmp->cur, target) && tmp->cur.pos != target.pos)
                {
                    nearst = &(tmp->cur);
                    dist = l2_distance(tmp->cur, target);
                }
                node* t;
                t = tmp->right;
                if (t != nullptr && is_in.find(t->cur.pos) == is_in.end())
                {
                    s.push(t);
                    is_in.insert(t->cur.pos);
                }
                t = tmp->left;
                if (t != nullptr && is_in.find(t->cur.pos) == is_in.end())
                {
                    s.push(t);
                    is_in.insert(t->cur.pos);
                }
            }
            else if (target.indices[cut_dim] < tmp->cur.indices[cut_dim])
            {
                node* t = tmp->left;
                if (t != nullptr && is_in.find(t->cur.pos) == is_in.end())
                {
                    s.push(t);
                    is_in.insert(t->cur.pos);
                }
            }
            else
            {
                node* t = tmp->right;
                if (t != nullptr && is_in.find(t->cur.pos) == is_in.end())
                {
                    s.push(t);
                    is_in.insert(t->cur.pos);
                }
            }
        }
        /*
                if (target.indices[cut_dim] <= tmp->cur.indices[cut_dim])
                    t = tmp->right;
                else
                    t = tmp->left;
                if (t != nullptr)
                    s.push(t);
            }
        }
        */
    }
    return nearst;

}
vector<point> kd_k_search(node* curnode, point target, int maxnum)
{

    float mu = 4.0f;
    point near = *(kd_nearst_search(curnode, target));
    float dist = l2_distance(near, target) * mu;
    stack<node*> s;
    set<int> is_in;
    vector<point> ans;
    priority_queue<pair<float, point>, vector<pair<float, point> >, pqcmp >q;
    while (curnode != nullptr)
    {
        s.push(curnode);
        is_in.insert(curnode->cur.pos);
        if (target.indices[curnode->depth % 3] <= (curnode->cur).indices[curnode->depth % 3])
            curnode = curnode->left;
        else
            curnode = curnode->right;
    }
    point* nearst = new point();
    /*
    point* nearst = &(s.top())->cur;
    ans.push_back(nearst);
    s.pop();
    */
    while (!s.empty())
    {
        node* tmp = s.top(); s.pop();
        if (tmp->left == nullptr && tmp->right == nullptr)
        {
            if (dist >= l2_distance(tmp->cur, target))
            {
                nearst = &(tmp->cur);
                if (nearst->pos == target.pos)
                {
                    ans.push_back(*nearst); continue;
                }
                if (q.size() < (unsigned int)maxnum)
                {

                    q.push(make_pair(l2_distance(*nearst, target), *nearst));
                }
                else
                {

                    q.pop();
                    q.push(make_pair(l2_distance(*nearst, target), *nearst));
                }
            }
        }
        else
        {
            int cut_dim = tmp->depth % 3;
            if (fabs(tmp->cur.indices[cut_dim] - target.indices[cut_dim]) < dist)
            {
                if (dist >= l2_distance(tmp->cur, target))
                {
                    nearst = &(tmp->cur);
                    if (q.size() < (unsigned int)maxnum)
                    {

                        q.push(make_pair(l2_distance(*nearst, target), *nearst));
                    }
                    else
                    {

                        q.pop();
                        q.push(make_pair(l2_distance(*nearst, target), *nearst));
                    }
                }
                node* t;
                /*
                if (target.indices[cut_dim] <= tmp->cur.indices[cut_dim])
                    t = tmp->right;
                else
                    t = tmp->left;
                if (t != nullptr )
                {
                    s.push(t);
                }
                */
                t = tmp->right;
                if (t != nullptr && is_in.find(t->cur.pos) == is_in.end())
                {
                    s.push(t);
                    is_in.insert(t->cur.pos);
                }
                t = tmp->left;
                if (t != nullptr && is_in.find(t->cur.pos) == is_in.end())
                {
                    s.push(t);
                    is_in.insert(t->cur.pos);
                }
            }
            else if (target.indices[cut_dim] < tmp->cur.indices[cut_dim])
            {
                node* t = tmp->left;
                if (t != nullptr && is_in.find(t->cur.pos) == is_in.end())
                {
                    s.push(t);
                    is_in.insert(t->cur.pos);
                }
            }
            else
            {
                node* t = tmp->right;
                if (t != nullptr && is_in.find(t->cur.pos) == is_in.end())
                {
                    s.push(t);
                    is_in.insert(t->cur.pos);
                }
            }
        }
    }
    while (!q.empty())
    {
        ans.push_back(q.top().second);
        q.pop();
    }
    is_in.clear();
    return ans;

}
bool is_right_point(point a, point b, point c, vector<point>& fa)
{
    /*
    if (l->indices[0] > r->indices[0])
    {
        point* t = l;l = r; r = t;
    }
    if (l->indices[0] == r->indices[0] && l->indices[1] > r->indices[1])
    {
        point* t = l; l = r; r = t;
    }
    if (l->indices[0] == r->indices[0] && l->indices[1] == r->indices[1] && l->indices[2] > r->indices[2])
    {
        point* t = l; l = r; r = t;
    }
    */
    point t1 = c - b;
    point t2 = a - b;
    point t = point(-1, t1.indices[1] * t2.indices[2] - t1.indices[2] * t2.indices[1], t1.indices[2] * t2.indices[0] - t1.indices[0] * t2.indices[2], t1.indices[0] * t2.indices[1] - t1.indices[1] * t2.indices[0]);
    t.norm();
    point cur_fa = fa[b.pos];
    float m = 0;
    for (int i = 0; i < 3; i++)
        m += cur_fa.indices[i] * t.indices[i];
    if (m > 0)
        return true;
    else
        return false;

}
float p_distance(point a, line* b)
{
    point* l = b->l;
    point* r = b->r;
    if (l->indices[0] > r->indices[0])
    {
        r = b->l; l = b->r;
    }
    float k = (r->indices[1] - l->indices[1]) / (r->indices[0] - l->indices[0]);
    float m = r->indices[1] - k * r->indices[0];
    float y1 = k * a.indices[0] + m;
    return fabs(a.indices[1] - y1) / sqrt(k * k + 1);
}
float get_angle(line* a, line* b)
{
    float x1 = a->r->indices[0] - a->l->indices[0], x2 = b->r->indices[0] - b->l->indices[0], y1 = a->r->indices[1] - a->l->indices[1], y2 = b->r->indices[1] - b->l->indices[1];
    float z1 = a->r->indices[2] - a->l->indices[2], z2 = b->r->indices[2] - b->l->indices[2];
    return acos((x1 * x2 + y1 * y2 + z1 * z2) / (sqrt(x1 * x1 + y1 * y1 + z1 * z1) * sqrt(x2 * x2 + y2 * y2 + z2 * z2)));
}
float get_angle_2d(line* a, line* b)
{
    float x1 = a->r->indices[0] - a->l->indices[0], x2 = b->r->indices[0] - b->l->indices[0], y1 = a->r->indices[1] - a->l->indices[1], y2 = b->r->indices[1] - b->l->indices[1];
    return acos((x1 * x2 + y1 * y2 ) / (sqrt(x1 * x1 + y1 * y1 ) * sqrt(x2 * x2 + y2 * y2 )+1e-8));
}
float get_angle(point* a, point* b, point* c)
{

    line* l1 = new line(a, b);
    line* l2 = new line(a, c);
    float angle1 = get_angle(l1, l2);
    if (fabs(angle1 - pi) < 1e-6)
        return (float)pi;
    *l1 = line(b, a); *l2 = line(b, c);
    float angle2 = get_angle(l1, l2);
    *l1 = line(c, a); *l2 = line(c, b);
    float angle3 = get_angle(l1, l2);
    delete l1;
    delete l2;
    return max(max(angle1, angle2), angle3);
}
float get_min_angle(point* a, point* b, point* c)
{
    line* l1 = new line(a, b);
    line* l2 = new line(a, c);
    float angle1 = get_angle(l1, l2);
    if (fabs(angle1 - pi) < 1e-6)
        return (float)pi;
    *l1 = line(b, a); *l2 = line(b, c);
    float angle2 = get_angle(l1, l2);
    *l1 = line(c, a); *l2 = line(c, b);
    float angle3 = get_angle(l1, l2);
    delete l1;
    delete l2;
    return min(min(angle1, angle2), angle3);
}
bool blur_triangle(triangle* t, float mmin, float mmax)
{
    float mi = get_min_angle(&t->points[0], &t->points[1], &t->points[2]), ma = get_angle(&t->points[0], &t->points[1], &t->points[2]);
    if (mi >= mmin && ma <= mmax)
        return true;
    return false;
}
bool is_in_triangle(triangle* tr, point* pt)
{

    line* l1 = new line(pt, &tr->points[0]);
    line* l2 = new line(pt, &tr->points[1]);
    line* l3 = new line(pt, &tr->points[2]);
    float sum = get_angle(l1, l2) + get_angle(l1, l3) + get_angle(l2, l3);
    delete l1;
    delete l2;
    delete l3;
    if (fabs(sum - 2 * pi) < 1e-6)
        return true;
    return false;

}
bool is_on_triangle(triangle* tr, point* pt)
{

    line* l[3];
    bool flag = false;
    int t = -1;
    for (int i = 0; i < 3; i++)
        l[i] = new line(pt, &tr->points[i]);
    for (int i = 0; i < 3; i++)
    {
        if (flag)
            break;
        if (fabs(get_angle(l[i], l[(i + 1) % 3]) - pi) < 1e-6)
            flag = true, t = i;
    }
    for (int i = 0; i < 3; i++)
        delete l[i];
    if (!flag)
        return false;
    float ll = tr->points[t].indices[0], r = tr->points[(t + 1) % 3].indices[0];
    if (pt->indices[0] > min(ll, r) && pt->indices[0] < max(ll, r))
        return true;
    return false;

}

bool cmp_ld(const point& a, const point& b, point ld_point)
{

    point* tmp_a = new point();
    point* tmp_b = new point();
    for (int i = 0; i < 3; i++)
    {
        tmp_a->indices[i] = a.indices[i] - ld_point.indices[i];
        tmp_b->indices[i] = b.indices[i] - ld_point.indices[i];
    }
    float x = tmp_a->indices[0] * tmp_b->indices[1] - tmp_a->indices[1] * tmp_b->indices[0];
    delete tmp_a;
    delete tmp_b;
    if (x > 0)
        return true;
    if (fabs(x) < 1e-12)
    {
        float s1 = 0, s2 = 0;
        for (int i = 0; i < 3; i++)
        {
            s1 += (a.indices[i] - ld_point.indices[i]) * (a.indices[i] - ld_point.indices[i]);
            s2 += (b.indices[i] - ld_point.indices[i]) * (b.indices[i] - ld_point.indices[i]);
        }
        if (fabs(s1 - s2) < 1e-12)
            return false;
        if (s1 < s2)
            return true;
    }
    return false;


}
point get_fa(vector<point>& points)
{

    int len = points.size();
    point* tmp_points = (point*)malloc(sizeof(point) * len);
    point* ans = new point();
    point* tmp = new point();
    for (int i = 0; i < len; i++)
    {
        for (int j = 0; j < 3; j++)
            (*tmp).indices[j] += points[i].indices[j];
    }
    for (int j = 0; j < 3; j++)
        (*tmp).indices[j] /= len;
    for (int i = 0; i < len; i++)
    {
        for (int j = 0; j < 3; j++)
            tmp_points[i].indices[j] = points[i].indices[j] - (*tmp).indices[j];
    }
    float tmp_mat[3][3] = { 0 };
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < len; k++)
                tmp_mat[i][j] += tmp_points[k].indices[i] * tmp_points[k].indices[j];
    free(tmp_points);
    Eigen::Matrix3f ma;
    ma << tmp_mat[0][0], tmp_mat[0][1], tmp_mat[0][2], tmp_mat[1][0], tmp_mat[1][1], tmp_mat[1][2], tmp_mat[2][0], tmp_mat[2][1], tmp_mat[2][2];
    Eigen::EigenSolver<Eigen::Matrix3f> es(ma);
    Eigen::Matrix3f D = es.pseudoEigenvalueMatrix();
    Eigen::Matrix3f V = es.pseudoEigenvectors();
    int j;
    float mim = 1e7;
    for (int i = 0; i < 3; i++)
        if (D(i, i) < mim)
        {
            mim = D(i, i);
            j = i;
        }
    for (int k = 0; k < 3; k++)
        ans->indices[k] = V(k, j);
    return *ans;

}
float point_to_plain(point* cur, point* fa, float bias)
{
    float t = 0;
    for (int i = 0; i < 3; i++)
        t += cur->indices[i] * fa->indices[i];
    t += bias;
    float div = 0;
    for (int i = 0; i < 3; i++)
        div += fa->indices[i] * fa->indices[i];
    return fabs(t) / sqrt(div);
}
bool is_in_circle(point p, triangle t)
{
    point tmp_p[2];
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 3; j++)
            tmp_p[i].indices[j] = t.points[i + 1].indices[j] - t.points[0].indices[j];
    float S = fabs(x_mul(tmp_p[0], tmp_p[1])) / 2;
    float mul = 1;
    for (int i = 0; i < 3; i++)
        mul *= l2_distance(t.points[i], t.points[(i + 1) % 3]);
    float r = mul / 4 / S;
    float xa = t.points[0].indices[0], ya = t.points[0].indices[1];
    float xb = t.points[1].indices[0], yb = t.points[1].indices[1];
    float xc = t.points[2].indices[0], yc = t.points[2].indices[1];
    float px[2], py[2], k[2], b[2];
    bool f[2] = { 0 };
    for (int i = 0; i < 2; i++)
    {
        px[i] = (t.points[0].indices[0] + t.points[i + 1].indices[0]) / 2;
        py[i] = (t.points[0].indices[1] + t.points[i + 1].indices[1]) / 2;
        if (t.points[i + 1].indices[0] == t.points[0].indices[0])
        {
            f[i] = true; b[i] = py[i]; continue;
        }
        k[i] = (t.points[0].indices[0] - t.points[1].indices[0]) / (t.points[i + 1].indices[1] - t.points[0].indices[1]);
        b[i] = py[i] - k[i] * px[i];
    }
    float x = (b[1] - b[0]) / (k[0] - k[1]);

    float y = k[0] * x + b[0];
    for (int i = 0; i < 2; i++)
        if (f[i])
        {
            x = b[i];
            y = k[i ^ 1] * x + b[i ^ 1];
        }
    point tmp = point(-1, x, y, 0);
    if (l2_distance(tmp, p) <= r)
        return true;
    return false;
}
void get_base(point fa, point* x_base, point* y_base)
{
    *x_base = point(-1, 1.0f, 1.0f, -(fa.indices[0] + fa.indices[1]) / fa.indices[2]);
    float y = fa.indices[0] * fa.indices[0] + fa.indices[2] * fa.indices[2] + fa.indices[0] * fa.indices[1];
    y /= (fa.indices[0] - fa.indices[1]);
    float x = (-fa.indices[2] * fa.indices[2] - fa.indices[1] * y) / fa.indices[0];
    *y_base = point(-1, x, y, fa.indices[2]);
    (*x_base).norm();
    (*y_base).norm();
}
point get_plain_point(point plain_fa, float d, point p, point fa)
{
    point ans = point();
    if (fa.indices[2] == 0)
    {
        ans.indices[2] = p.indices[2];
        if (fa.indices[1] == 0)
        {
            ans.indices[1] = p.indices[1];
            ans.indices[0] = (-d - plain_fa.indices[1] * ans.indices[1] - plain_fa.indices[2] * ans.indices[2]) / plain_fa.indices[0];
            return ans;
        }
        else
        {
            ans.indices[1] = (-d - plain_fa.indices[1] * ans.indices[1] - plain_fa.indices[2] * ans.indices[2] + plain_fa.indices[0] * fa.indices[0] / fa.indices[1] * p.indices[1]) / (plain_fa.indices[0] * fa.indices[0] / fa.indices[1] + plain_fa.indices[1]);
            ans.indices[0] = fa.indices[0] / fa.indices[1] * (ans.indices[1] - p.indices[1]) + p.indices[0];
            return ans;
        }
    }
    float tmp = plain_fa.indices[0] * (fa.indices[0] / fa.indices[2]) * p.indices[2] - plain_fa.indices[0] * p.indices[0] - plain_fa.indices[1] * p.indices[1] + plain_fa.indices[1] * fa.indices[1] / fa.indices[2] * p.indices[2] - d;
    ans.indices[2] = tmp / (plain_fa.indices[0] * fa.indices[0] / fa.indices[2] + plain_fa.indices[1] * fa.indices[1] / fa.indices[2] + plain_fa.indices[2]);
    ans.indices[1] = p.indices[1] + fa.indices[1] / fa.indices[2] * (p.indices[2] - fa.indices[2]);
    ans.indices[0] = p.indices[0] + fa.indices[0] / fa.indices[2] * (p.indices[2] - fa.indices[2]);
    return ans;
}
float dot_mul(point a, point b)
{
    float ans = 0;
    for (int i = 0; i < 3; i++)
        ans += a.indices[i] * b.indices[i];
    return ans;
}
bool cmp_zero(point a, point b, point zero)
{
    point ta = a - zero, tb = b - zero;
    float l1 = l2_distance(a, zero), l2 = l2_distance(b, zero);
    ta.norm(); tb.norm();
    //return (ta.indices[1] * ta.indices[1]) < (tb.indices[1] * tb.indices[1]);
    float angle_a = atan(ta.indices[1] / ta.indices[0]), angle_b = atan(tb.indices[1] / tb.indices[0]);
    if (ta.indices[0] < 0)
        angle_a = (float)pi + angle_a;
    if (angle_a < 0)
        angle_a = (float)pi * 2 + angle_a;
    if (tb.indices[0] < 0)
        angle_b = (float)pi + angle_b;
    if (angle_b < 0)
        angle_b = (float)pi * 2 + angle_b;
    if (angle_a == angle_b)
        return l1 < l2;
    return angle_a < angle_b;

}
int check_if_triangle(point* points, int pre, int i, int late, int terminal, int len)
{
    if (i == terminal)
        return 1;
    if (pre == -1)
        pre = len - 1;
    point mid_a = points[pre] / 2, mid_b = points[i] / 2, mid_c = points[late] / 2;
    float k_a = points[pre].indices[1] / points[pre].indices[0];
    float k_b = points[i].indices[1] / points[i].indices[0];
    float k_c = points[late].indices[1] / points[late].indices[0];
    float x_ac = (k_a * mid_c.indices[0] + (mid_c.indices[1] - mid_a.indices[1]) * k_a * k_c - k_c * mid_a.indices[0]) / (k_a - k_c);
    float y_ac = -1 / k_a * (x_ac - mid_a.indices[0]) + mid_a.indices[1];
    float k = -1 / k_b;
    float b = y_ac - k * x_ac;
    bool flag = true;
    float y = k * mid_b.indices[0] + b;
    if (b < 0 && y > mid_b.indices[1])
        flag = false;
    if (b > 0 && y < mid_b.indices[1])
        flag = false;
    if (flag)
        return 1;
    else
    {
        check_if_triangle(points, pre - 1, pre, late, terminal, len);
        return 0;
    }
}
int find_pa(vector<int>& pa, int pos)
{
    /*
    if (pa[pos] > int(pa.size()))
    {
        int t = pos;

    }
    */
    return (pa[pos] == pos) ? pos : pa[pos] = find_pa(pa, pa[pos]);
}
void smallest_tree(vector<point>& fa, point* points, node* root)
{
    int l = fa.size();
    int* from = (int*)malloc(sizeof(int) * l * 20);
    int* tov = (int*)malloc(sizeof(int) * l * 20);
    memset(from, -1, sizeof(int) * l * 20);
    memset(tov, -1, sizeof(int) * l * 20);
    int tot = 0;
    priority_queue<pq_node, vector<pq_node>, greater<pq_node> > q;
    for (int i = 0; i < l; i++)
    {
        vector<point>tmp = kd_k_search(root, points[i], 20);
        int ll = tmp.size();
        for (int j = 0; j < ll; j++)
        {
            if (tmp[j].pos == i)
                continue;
            from[tot] = i; tov[tot] = tmp[j].pos;
            q.push(pq_node(1 - fabs(dot_mul(fa[i], fa[tmp[j].pos])), tot)); tot++;
        }
    }
    tot = 0;
    int* head = (int*)malloc(sizeof(int) * l);
    memset(head, 0, sizeof(int) * l);
    int* too = (int*)malloc(sizeof(int) * l);
    int* nxt = (int*)malloc(sizeof(int) * l);
    vector<int>pa;
    for (int i = 0; i < l; i++)
        pa.push_back(i);
    while (!q.empty())
    {
        pq_node t = q.top(); q.pop();
        int u = from[t.pos], v = tov[t.pos];
        if (u > v)
        {
            int tmp = u; u = v; v = tmp;
        }
        int pu = find_pa(pa, u), pv = find_pa(pa, v);
        if (pu != pv)
        {
            pa[pv] = pu;
            nxt[++tot] = head[pu];
            head[pu] = tot;
            too[tot] = pv;
        }
    }
    int pq_root = pa[0];
    //cout << pa[0] << endl;
    queue<int> qu;
    qu.push(pq_root);
    fa[pq_root] = point(0, 0, 0, 0) - fa[pq_root];
    while (!qu.empty())
    {
        int t = qu.front(); qu.pop();
        for (int j = head[t]; j != 0; j = nxt[j])
        {
            int v = too[j];
            if (dot_mul(fa[t], fa[v]) < 0)
            {
                fa[v] = point(0, 0, 0, 0) - fa[v];
            }
            qu.push(v);
        }
    }
    free(from); free(tov); free(head); free(too); free(nxt);
}