
#include "greedy_triangle.h"
using namespace std;



list<triangle> triangle_cut(vector<point>& points, vector<point>& fa, point near, map<int, vector<point> >& ma_v)
{
    int len = points.size();
    list<triangle> ans;
    /*
    if (len <= 2)
        return ans;
        */
    list<point> othpoint;
    list<point> hullpoint;
    point cur_fa = fa[points[0].pos];
    cur_fa.norm();

    float mam = 0;
    for (int i = 0; i < 3; i++)
        mam = max(mam, fabs((near - points[0]).indices[i]));
    point x_base, y_base;
    get_base(cur_fa, &x_base, &y_base);
    float bias = 0;

    int mq_size = min(20, (int)points.size() - 1);
    priority_queue<pair<float, point>, vector<pair<float, point> >, pqcmp >q;
    for (int i = 1; i < len; i++)
    {
        float t = l2_distance(points[i], points[0]);
        if ((int)q.size() >= mq_size)
        {
            float ma = q.top().first;
            if (t < ma)
            {
                q.pop();
                q.push(make_pair(t, points[i]));
            }
        }
        else
            q.push(make_pair(t, points[i]));
    }
    point mq_point[20];
    int mq_tot = 0;
    while (!q.empty())
    {
        mq_point[mq_tot++] = q.top().second; q.pop();
    }
    point v1 = (near - points[0]) - fa[near.pos] * dot_mul(cur_fa, (near - points[0]));
    point vi[20];
    float kpvi[20];
    pq_node angles[20];
    int mq_pos = -1;
    point mq_ver;
    for (int i = 0; i < mq_size; i++)
    {
        point tmp = (mq_point[i] - points[0]);
        vi[i] = tmp - fa[mq_point[i].pos] * dot_mul(cur_fa, tmp);
        kpvi[i] = 2 * dot_mul(cur_fa - fa[mq_point[i].pos], tmp) / (tmp.length() * tmp.length());
        angles[i].dis = acos(dot_mul(vi[i], v1) / (vi[i].length() * v1.length()));
        angles[i].pos = i;
        if (mq_pos == -1 && angles[i].dis < pi / 2)
        {
            mq_pos = i;
            float check1 = v1.length() * v1.length();
            float check2 = dot_mul(vi[i], v1) / check1;
            mq_ver = vi[i] - v1 * check2;
            //mq_ver = vi[i] - v1 * dot_mul(vi[i], v1) / (v1.length() * v1.length());
            mq_ver.norm();
        }
    }
    for (int i = 0; i < mq_size; i++)
    {
        if (dot_mul(vi[i], mq_ver) < 0)
            angles[i].dis = 2 * (float)pi - angles[i].dis;
    }
    sort(angles, angles + mq_size);
    float mq[3][3] = { 0 };
    for (int i = 0; i < mq_size; i++)
    {
        float tmp = angles[(i + 1) % 20].dis - angles[(i + 20) % 20].dis;
        if (tmp < 0)
            tmp += 2 * float(pi);
        tmp = tmp / 4 / (float)pi * kpvi[angles[i].pos];
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
            {
                mq[j][k] += vi[angles[i].pos].indices[j] * vi[angles[i].pos].indices[k];
            }
    }
    Eigen::Matrix3f matrix;
    matrix << mq[0][0], mq[0][1], mq[0][2], mq[1][0], mq[1][1], mq[1][2], mq[2][0], mq[2][1], mq[2][2];
    Eigen::EigenSolver<Eigen::Matrix3f> es(matrix);
    Eigen::Matrix3f D = es.pseudoEigenvalueMatrix();
    Eigen::Matrix3f V = es.pseudoEigenvectors();
    point pc[3];
    float kc[3];
    int mi;
    float mim = 1e7;
    for (int i = 0; i < 3; i++)
    {
        if (D(i, i) < mim)
        {
            mim = D(i, i);
            mi = i;
        }
    }

    for (int tot = 0, i = 0; i < 3; i++)
    {
        if (i != mi)
        {
            kc[tot] = D(i, i);
            for (int j = 0; j < 3; j++)
            {
                pc[tot].indices[j] = V(j, i);
            }
            tot++;
        }
    }
    float ktmp = kc[0]; kc[0] = 3 * kc[1] - kc[0]; kc[1] = 3 * ktmp - kc[1];
    point ptmp = pc[0]; pc[0] = pc[1] * 3 - pc[0]; pc[1] = ptmp * 3 - pc[1];
    if (kc[0] > kc[1])
    {
        ktmp = kc[0]; kc[0] = kc[1]; kc[1] = ktmp;
        ptmp = pc[0]; pc[0] = pc[1]; pc[1] = ptmp;
    }
    pc[0].norm(); pc[1].norm();
    float mp = 2 * kc[1] / kc[0], near_dis = l2_distance(near, points[0]);
    float delta = 2 * dot_mul(cur_fa - fa[near.pos], (near - points[0])) / ((near - points[0]).length());
    float delta_thres = (sqrt(1 + 4 * delta * delta) - 1) / kc[0];

    vector<point> ma_vp;
    ma_vp.push_back(pc[0]);
    ma_vp.push_back(pc[1]);
    ma_v[points[0].pos] = ma_vp;

    for (int i = 0; i < 3; i++)
        bias -= points[0].indices[i] * cur_fa.indices[i];
    float threshold = l2_distance(near, (points[0])) * 1.5f;
    float d = 0;
    for (int i = 0; i < 3; i++)
        d += cur_fa.indices[i] * points[0].indices[i];
    point v = near - points[0];
    point tmp_fa = cur_fa - fa[near.pos];
    float thres = dot_mul(tmp_fa, v) / (v.length() + (float)1e-9);

    mp = fabs(mp);
    //cout << dot_mul(cur_fa, pc[0]) << ' ' << dot_mul(cur_fa, pc[1]) << endl;
    for (int i = 0; i < len; i++)
    {

        if (l2_distance(points[i], points[0]) > mp * near_dis)
            continue;
        float x = dot_mul(points[i] - points[0], pc[0]);
        float y = dot_mul(points[i] - points[0], pc[1]);
        float z = dot_mul(points[i] - points[0], cur_fa);
        if (fabs(z) > delta_thres)
            continue;

        //if (x < 0)continue;

        point tmp_p = points[i] - points[0];
        float t_ma = 0;
        for (int j = 0; j < 3; j++)
            t_ma = max(t_ma, fabs(tmp_p.indices[j]));
        if (t_ma > 2 * mam)
            continue;


        point* tmp = new point(points[i].pos, x, y, 0.0f);


        //point plain_point = get_plain_point(cur_fa,d,points[i], fa[points[i].pos]);
        //point miu = plain_point - points[0];
        //point* tmp = new point(points[i].pos,dot_mul(miu, x_base), dot_mul(miu, y_base), 0.0f);

        //float t = -points[i]->indices[2] / fa[points[i]->pos].indices[2];
        //point* tmp = new point(points[i]->pos, t * fa[points[i]->pos].indices[0] + points[i]->indices[0], t * fa[points[i]->pos].indices[1] + points[i]->indices[1], 0);
        othpoint.push_back(*tmp);
        delete tmp;
    }

    len = othpoint.size();
    if (len <= 2)
        return ans;
    if (len == 3)
    {
        triangle tmp;
        int i = 0;
        for (list<point>::iterator it = othpoint.begin(); it != othpoint.end(); i++, it++)
            tmp.points[i] = (*it);
        ans.push_back(tmp);
        return ans;
    }

    point ld_point = othpoint.front();
    point* tmp_points = (point*)malloc(sizeof(point) * len);
    int tmp_points_i = 0;
    othpoint.pop_front();
    for (list<point>::iterator m = othpoint.begin(); m != othpoint.end(); m++)
        tmp_points[tmp_points_i++] = *m;
    for (int i = 0; i < tmp_points_i; i++)
        for (int j = i + 1; j < tmp_points_i; j++)
        {
            if (!cmp_zero(tmp_points[i], tmp_points[j], ld_point))
            {
                point t = tmp_points[i]; tmp_points[i] = tmp_points[j]; tmp_points[j] = t;
            }
        }
    othpoint.clear();
    map<int, int> ma;
    int pflag = near.pos;
    int near_pos = 0;
    for (; near_pos < tmp_points_i; near_pos++)
        if (tmp_points[near_pos].pos == pflag)
            break;

    if (near_pos == tmp_points_i)
    {
        near = point(0, 1, 1, 1);
        near_pos = 0;
        for (int i = 0; i < tmp_points_i; i++)
            if (l2_distance(tmp_points[i], ld_point) < l2_distance(near, ld_point))
                near = tmp_points[i], near_pos = i;
    }

    int i = near_pos + 1;
    ma[near_pos] = 1;

    while (i != near_pos)
    {

        if (i == tmp_points_i)
            i = 0;
        int pre_i = i - 1;
        if (pre_i == -1)
            pre_i = tmp_points_i - 1;
        int late_i = i + 1;
        if (late_i == tmp_points_i)
            late_i = 0;
        if (i == near_pos)
            break;
        ma[i] = check_if_triangle(tmp_points, pre_i, i, late_i, near_pos, tmp_points_i);
        i++;
    }
    for (int i = 0; i < tmp_points_i; i++)
        if (ma[i] == 1)
            othpoint.push_back(tmp_points[i]);

    for (list<point>::iterator it = othpoint.begin(); it != othpoint.end(); it++)
    {
        list<point>::iterator nxt_it = it;
        nxt_it++;
        /*
        if (nxt_it == othpoint.end())
            nxt_it = othpoint.begin();
         */
        if (nxt_it == othpoint.end())
            break;
        triangle tmp_ans = triangle(ld_point, *it, *nxt_it);

        if (!blur_triangle(&tmp_ans, (float)pi / 10, (float)pi))
            continue;

        ans.push_back(tmp_ans);
    }


    othpoint.clear();
    hullpoint.clear();
    return ans;

}


void check_hole(int* is_hole, vector<face>& ans, vector<point>& cur_point, vector<point>& fa, map<pair<int, int>, int>& ma_line, set<int>& is_in_hole)
{
    int l = cur_point.size();
    //cout << l << endl;
    if (l <= 2 || l >= 20)
        return;
    for (int i = 0; i < l; i++)
        if (is_hole[cur_point[i].pos] == 0)
            is_hole[cur_point[i].pos] = 1;
    point cur_fa = point();
    for (int i = 0; i < l; i++)
        cur_fa = cur_fa + fa[cur_point[i].pos];
    cur_fa = cur_fa / l;
    for (int i = 0; i < l; i++)
    {
        //cout << ans.size() << endl;
        /*
        if (is_in_hole.find(cur_point[i].pos) == is_in_hole.end())
        {
            is_in_hole.insert(cur_point[i].pos);
        }
        */
        //ans.push_back(face(cur_point[i].pos, cur_point[(i + 1) % l].pos, cur_point[((i + 1) % l + 1) % l].pos));
    }

}
void dfs(int u, int* head, int* nxt, int* to, int* vis, int* is_hole, vector<int>& hole_s, vector<face>& ans, point* points, vector<point>& fa, map<pair<int, int>, int>& ma_line, set<int>& is_in_hole)
{
    vis[u] = 1;
    hole_s.push_back(u);
    for (int k = head[u]; k; k = nxt[k])
    {
        int v = to[k];
        if (vis[v] == 1)
        {
            vector<point> cur_point;
            int len = hole_s.size();
            for (int j = len - 1; j >= 0; j--)
            {
                cur_point.push_back(points[hole_s[j]]);
                if (hole_s[j] == v)
                {
                    break;
                }
            }
            check_hole(is_hole, ans, cur_point, fa, ma_line, is_in_hole);
            continue;
        }
        if (vis[v] == 0)
            dfs(v, head, nxt, to, vis, is_hole, hole_s, ans, points, fa, ma_line, is_in_hole);
    }

    /*
    if (is_hole[u] == 0)
    {
        vector<point> cur_point;
        int len = hole_s.size();
        for (int j = len - 1; j >= 0; j--)
        {
            cur_point.push_back(points[hole_s[j]]);
        }
        cout << cur_point.size() << endl;
        check_hole(is_hole,ans, cur_point, fa);
    }
    */
    vis[u] = 2;
    hole_s.pop_back();
}
bool check_line(point* points, vector<point>& fa, pair<int, int> tmp_line, pair<int, int> tmp_check_line, map<int, vector<point> >& ma_v)
{
    int u = tmp_line.first, v = tmp_line.second;
    int u1 = tmp_check_line.first, v1 = tmp_check_line.second;
    point cur_fa = fa[u];
    cur_fa = cur_fa.norm();
    cur_fa.norm();
    if (ma_v.find(u) == ma_v.end())
        return false;
    point x_base = ma_v[u][0];
    x_base.norm();
    point y_base = ma_v[u][1];
    y_base.norm();
    point u_cap = point(points[u].pos, 0, 0, 0);
    point v_cap = point(points[v].pos, dot_mul(x_base, points[v] - points[u]), dot_mul(y_base, points[v] - points[u]), dot_mul(cur_fa, points[v] - points[u]));
    point u1_cap = point(points[u1].pos, dot_mul(x_base, points[u1] - points[u]), dot_mul(y_base, points[u1] - points[u]), dot_mul(cur_fa, points[u1] - points[u]));
    point v1_cap = point(points[v1].pos, dot_mul(x_base, points[v1] - points[u]), dot_mul(y_base, points[v1] - points[u]), dot_mul(cur_fa, points[v1] - points[u]));
    //cout << v_cap.indices[2] << ' ' << u1_cap.indices[2] << ' ' << v1_cap.indices[2] << endl;
    float c1 = x_mul(v_cap - u_cap, u1_cap - u_cap), c2 = x_mul(v_cap - u_cap, v1_cap - u_cap);
    float c3 = x_mul(u_cap - u1_cap, v1_cap - u1_cap), c4 = x_mul(v1_cap - u1_cap, v_cap - u1_cap);
    if (c1 * c2 < 0 && c3 * c4 < 0)
        return true;
    return false;
}
bool check_line_inter(pair<int, int> cur_line, pair<int, int> pre_line, point* points)
{
    point cur = points[cur_line.second] - points[cur_line.first];
    point pre = points[pre_line.second] - points[pre_line.first];
    point f = x_mul3d(cur, pre);
    if (f.length() < 1e-9)
        return true;
    point t = points[cur_line.first] - points[pre_line.first];
    float is_plain = dot_mul(t, f);
    if (fabs(is_plain) < 1e-9)
        return true;
    float t2 = dot_mul(x_mul3d(t, cur), f) / (f.length() * f.length());
    point center = points[pre_line.first] + pre * t2;
    float c1 = dot_mul(points[pre_line.second] - center, points[pre_line.first] - center), c2 = dot_mul(points[cur_line.second] - center, points[cur_line.first] - center);
    if (c1 < 0 && c2 < 0)
        return false;
    return true;
}
void check_dfs_hole(vector<point>& ma_v, vector<point>& fa, set<face>& s, vector<face>& ans, vector<point>& cur_point, int target)
{
    int len = cur_point.size();
    if (len <= 2)
        return;
    int pos = -1;
    for (int i = 0; i < len; i++)
        if (cur_point[i].pos == target)
        {
            pos = i;
            break;
        }
    if (pos == -1)
        return;
    float mam = 0;

    point tmp_p = cur_point[pos]; cur_point[pos] = cur_point[0]; cur_point[0] = tmp_p;

    point cur_fa = fa[target];
    int near = -1;
    float minm = 1e7;
    for (int i = 0; i < len; i++)
    {
        if (i != pos && l2_distance(cur_point[i], cur_point[pos]) < minm)
        {
            near = i;
            minm = l2_distance(cur_point[i], cur_point[pos]);
        }
    }

    for (int i = 0; i < 3; i++)
        mam = max(mam, fabs((cur_point[near] - cur_point[0]).indices[i]));


    list<point> othpoint;
    point near_p = cur_point[near];
    for (int i = 0; i < len; i++)
    {
        point tmp_p = cur_point[i] - cur_point[0];
        float t_ma = 0;
        for (int j = 0; j < 3; j++)
            t_ma = max(t_ma, fabs(tmp_p.indices[j]));
        /*
        if (t_ma > 3 * mam)
            continue;
        */
        float x = dot_mul(cur_point[i] - cur_point[0], ma_v[0]);
        float y = dot_mul(cur_point[i] - cur_point[0], ma_v[1]);
        float z = dot_mul(cur_point[i] - cur_point[0], cur_fa);
        point* tmp = new point(cur_point[i].pos, x, y, 0.0f);
        othpoint.push_back(*tmp);
        delete tmp;
    }
    int l = othpoint.size();
    cout << l << endl;
    if (l <= 2)
        return;
    if (l == 3)
    {
        triangle tmp;
        int i = 0;
        for (list<point>::iterator it = othpoint.begin(); it != othpoint.end(); i++, it++)
            tmp.points[i] = (*it);
        face t = face(&tmp);
        if (s.find(t) == s.end())
        {
            ans.push_back(t);
            s.insert(t);
        }
        return;
    }
    point ld_point = othpoint.front();
    point* tmp_points = (point*)malloc(sizeof(point) * len);
    int tmp_points_i = 0;
    othpoint.pop_front();
    for (list<point>::iterator m = othpoint.begin(); m != othpoint.end(); m++)
        tmp_points[tmp_points_i++] = *m;
    for (int i = 0; i < tmp_points_i; i++)
        for (int j = i + 1; j < tmp_points_i; j++)
        {
            if (!cmp_zero(tmp_points[i], tmp_points[j], ld_point))
            {
                point t = tmp_points[i]; tmp_points[i] = tmp_points[j]; tmp_points[j] = t;
            }
        }
    othpoint.clear();
    map<int, int> ma;
    for (int i = 0; i < tmp_points_i; i++)
        if (tmp_points[i].pos == near_p.pos)
        {
            near = i;
            break;
        }
    int i = near + 1;
    ma[near] = 1;
    while (i != near)
    {

        if (i == tmp_points_i)
            i = 0;
        int pre_i = i - 1;
        if (pre_i == -1)
            pre_i = tmp_points_i - 1;
        int late_i = i + 1;
        if (late_i == tmp_points_i)
            late_i = 0;
        if (i == near)
            break;
        ma[i] = check_if_triangle(tmp_points, pre_i, i, late_i, near, tmp_points_i);
        i++;
    }

    for (int i = 0; i < tmp_points_i; i++)
        if (ma[i] == 1)
            othpoint.push_back(tmp_points[i]);


    for (list<point>::iterator it = othpoint.begin(); it != othpoint.end(); it++)
    {
        list<point>::iterator nxt_it = it;
        nxt_it++;
        /*
        if (nxt_it == othpoint.end())
            nxt_it = othpoint.begin();
         */
        if (nxt_it == othpoint.end())
            break;
        triangle tmp_ans = triangle(ld_point, *it, *nxt_it);
        /*
        if (!blur_triangle(&tmp_ans, pi / 18, pi ))
            continue;
        */
        face t = face(&tmp_ans);
        if (s.find(t) == s.end())
            s.insert(t), ans.push_back(t);
    }

}
void check_dfs(int u, point* points, vector<point>& ma_v, vector<point>& fa, set<face>& s, vector<face>& ans, int* head, int* nxt, int* to, int* is_vis, vector<int>& hole_s, int depth, int target)
{
    is_vis[u] = 1;
    if (depth >= 10)
    {
        is_vis[u] = 2;
        return;
    }
    hole_s.push_back(u);
    for (int k = head[u]; k; k = nxt[k])
    {
        int v = to[k];
        if (is_vis[v] == 1)
        {
            vector<point> cur_point;
            int len = hole_s.size();
            for (int j = len - 1; j >= 0; j--)
            {
                cur_point.push_back(points[hole_s[j]]);
                if (hole_s[j] == v)
                {
                    break;
                }
            }
            check_dfs_hole(ma_v, fa, s, ans, cur_point, target);
            continue;
        }
        if (is_vis[v] == 0)
            check_dfs(v, points, ma_v, fa, s, ans, head, nxt, to, is_vis, hole_s, depth + 1, target);
    }
    is_vis[u] = 2;
    hole_s.pop_back();

}
vector<face> greedy_triangle(vector<float>& ini, vector<point>& fa, vector<int>& point_color)
{

    int l = ini.size();
    set<face> s;
    vector<face> ans;
    point* points = (point*)malloc(sizeof(point) * (l / 3));
    map<pair<int, int>, int> ma_line;
    for (int i = 0; i < l / 3; i++)
    {
        points[i] = point(i, ini[i * 3], ini[i * 3 + 1], ini[i * 3 + 2]);
    }
    node* root = kd_build(points, l / 3, 0);
#pragma omp parallel for
    for (int i = 0; i < l / 3; i++)
    {
        //printf("get_fa %d\n", i);
        vector<point> tmp = kd_k_search(root, points[i], 10);
        point t = get_fa(tmp);
        t.pos = points[i].pos;
        t.norm();
        /*
        float check = 0;
        for (int j = 0; j < 3; j++)
            check += t.indices[j] * points[i].indices[j];
        if (check > 0)
            for (int i = 0; i < 3; i++)
                t.indices[i] = -t.indices[i];
        */
#pragma omp critical
        {
            fa.push_back(t);
        }

    }
    sort(fa.begin(), fa.end(), cmppos);
    smallest_tree(fa, points, root);

    map<int, int> near_vis;
    sort(points, points + l / 3, cmppos);
    map<pair<pair<int, int>, int>, vector<int> > ma_point;
    set<pair<int, int> >is_not_hole;
    map<int, vector<point> > ma_v;
    vector<pair<int, int> > vec_pa;
    map<pair<int, int>, int> ma_pa;
    node* edge_root = kd_build(&point(), 1, 0);

#pragma omp parallel for
    for (int i = 0; i < l / 3; i++)
    {

        //if (i != 18768)continue;
        //cout << i << endl;

        vector<point> tmp = kd_k_search(root, points[i], 10);
        list<triangle> tmp_triangle;
        /*
        _try {
            tmp_triangle = triangle_cut(tmp, *tmp_point);
        }
        _finally
        {
            cout << i << endl;
        }
        */
        //printf("%d %d\n", i,tmp.size());
        /*
        if(i == 22373)
            tmp_triangle = triangle_cut(tmp, fa, kd_nearst_search(root, points[i]));
         */
        for (int j = 0; j < (int)tmp.size(); j++)
        {
            if (tmp[j].pos == points[i].pos)
            {
                point tmp_p = tmp[0];
                tmp[0] = tmp[j];
                tmp[j] = tmp_p;
            }
        }
        tmp_triangle = triangle_cut(tmp, fa, *kd_nearst_search(root, points[i]), ma_v);
        bool flags[20];
        for (int j = 0; j < 20; j++)
            flags[j] = true;
        int j = 0;
        for (list<triangle>::iterator it = tmp_triangle.begin(); it != tmp_triangle.end(); it++, j++)
        {
            face t = face(&(*it));
            point p_center[3];
            for (int k = 0; k < 3; k++)
            {
                int u = min(t.indices[k], t.indices[(k + 1) % 3]), v = max(t.indices[k], t.indices[(k + 1) % 3]);
                pair<int, int> tmp_line = make_pair(u, v);
                p_center[k] = points[u] + points[v];
                p_center[k] = p_center[k] / 2;
                p_center[k].pos = -1;
                vector<point> tmp_p = kd_k_search(edge_root, p_center[k], 10);
                for (int m = 0; m < tmp_p.size(); m++)
                {
                    if (tmp_p[m].pos == -1)
                        continue;

                    flags[j] = check_line_inter(tmp_line, vec_pa[tmp_p[m].pos], points);
                }
            }
        }

#pragma omp critical
        {
            /*
            int v = (*kd_nearst_search(root, points[i])).pos;
            if (near_vis.find(v) != near_vis.end())
            {
                if (near_vis[v] == i)
                    continue;
            }
            near_vis[i] = v;
            */
            int j = 0;
            for (list<triangle>::iterator it = tmp_triangle.begin(); it != tmp_triangle.end(); it++, j++)
            {
                face t = face(&(*it));
                if (s.find(t) == s.end())
                {
                    bool flag = flags[j];
                    for (int k = 0; k < 3; k++)
                    {
                        int u = min(t.indices[k], t.indices[(k + 1) % 3]), v = max(t.indices[k], t.indices[(k + 1) % 3]);
                        pair<int, int> tmp_line = make_pair(u, v);

                        if (ma_line.find(tmp_line) == ma_line.end())
                        {
                            continue;
                        }
                        else
                        {
                            if (ma_line[tmp_line] == 2)
                            {
                                flag = false;
                                break;
                            }
                            if (ma_line[tmp_line] != int(is_right_point(points[t.indices[(k + 2) % 3]], points[u], points[v], fa)))
                            {
                                ma_line[tmp_line] = 2;
                            }
                            else
                            {
                                flag = false;
                                break;
                            }
                        }
                    }
                    if (flag)
                    {
                        if (!flag)
                            continue;
                        point p_center[3];
                        for (int k = 0; k < 3; k++)
                        {
                            int u = min(t.indices[k], t.indices[(k + 1) % 3]), v = max(t.indices[k], t.indices[(k + 1) % 3]);
                            pair<int, int> tmp_line = make_pair(u, v);
                            if (ma_pa.find(tmp_line) == ma_pa.end())
                            {
                                ma_pa[tmp_line] = vec_pa.size();
                                vec_pa.push_back(tmp_line);
                                p_center[k] = points[u] + points[v];
                                p_center[k] = p_center[k] / 2;
                                p_center[k].pos = ma_pa[tmp_line];

                                kd_add(edge_root, p_center[k]);
                            }
                        }

                        for (int k = 0; k < 3; k++)
                        {
                            int u = min(t.indices[k], t.indices[(k + 1) % 3]), v = max(t.indices[k], t.indices[(k + 1) % 3]);
                            pair<int, int> tmp_line = make_pair(u, v);
                            if (ma_line.find(tmp_line) == ma_line.end())
                            {
                                ma_line[tmp_line] = is_right_point(points[t.indices[(k + 2) % 3]], points[u], points[v], fa);
                            }
                        }


                        ans.push_back(t); s.insert(t);
                    }
                }
            }
        }
    }

    ma_line.clear();
    vector<pair<int, int> > hole_line;


    for (int i = 0; i < int(ans.size()); i++)
    {
        face t = ans[i];
        for (int k = 0; k < 3; k++)
        {
            int u = min(t.indices[k], t.indices[(k + 1) % 3]), v = max(t.indices[k], t.indices[(k + 1) % 3]);
            pair<int, int> tmp_line = make_pair(u, v);
            if (ma_line.find(tmp_line) == ma_line.end())
            {
                ma_line[tmp_line] = is_right_point(points[t.indices[(k + 2) % 3]], points[u], points[v], fa);
                vector<int> tmp_vector;
                vector<int> tmp_vector_1;
                tmp_vector.push_back(t.indices[(k + 2) % 3]);
                ma_point[make_pair(tmp_line, ma_line[tmp_line])] = tmp_vector;
                ma_point[make_pair(tmp_line, !ma_line[tmp_line])] = tmp_vector_1;
            }
            else
            {
                if (ma_line[tmp_line] != int(is_right_point(points[t.indices[(k + 2) % 3]], points[u], points[v], fa)))
                {
                    if (is_not_hole.find(tmp_line) == is_not_hole.end())
                        is_not_hole.insert(tmp_line);
                    ma_point[make_pair(tmp_line, !ma_line[tmp_line])].push_back(t.indices[(k + 2) % 3]);
                }
                else
                {
                    ma_point[make_pair(tmp_line, ma_line[tmp_line])].push_back(t.indices[(k + 2) % 3]);
                }
            }
        }
    }
    for (int i = 0; i < int(ans.size()); i++)
    {
        face t = ans[i];
        //cout << i << endl;
        for (int k = 0; k < 3; k++)
        {
            int u = min(t.indices[k], t.indices[(k + 1) % 3]), v = max(t.indices[k], t.indices[(k + 1) % 3]);
            pair<int, int> tmp_line = make_pair(u, v);
            if (is_not_hole.find(tmp_line) != is_not_hole.end())
            {
                vector<int> a1 = ma_point[make_pair(tmp_line, 0)];
                vector<int> a2 = ma_point[make_pair(tmp_line, 1)];
                int l1 = a1.size(), l2 = a2.size();
                for (int x = 0; x < l1; x++)
                    for (int y = 0; y < l2; y++)
                    {
                        int u1 = min(a1[x], a2[y]), v1 = max(a1[x], a2[y]);
                        pair<int, int> tmp_check_line = make_pair(u1, v1);
                        bool flag = check_line(points, fa, tmp_line, tmp_check_line, ma_v);
                        if (flag)
                        {
                            if (is_not_hole.find(tmp_check_line) == is_not_hole.end())
                                is_not_hole.insert(tmp_check_line);
                        }
                    }
            }
        }
    }
    ma_point.clear();
    for (int i = 0; i < int(ans.size()); i++)
    {
        face t = ans[i];
        for (int k = 0; k < 3; k++)
        {
            int u = min(t.indices[k], t.indices[(k + 1) % 3]), v = max(t.indices[k], t.indices[(k + 1) % 3]);
            pair<int, int> tmp_line = make_pair(u, v);
            if (is_not_hole.find(tmp_line) == is_not_hole.end())
                hole_line.push_back(tmp_line);
        }
    }

    int* head = (int*)malloc(sizeof(int) * l);
    memset(head, 0, sizeof(int) * l);
    int len_line = hole_line.size();
    int* nxt = (int*)malloc(sizeof(int) * (len_line + 1) * 2);
    int* to = (int*)malloc(sizeof(int) * (len_line + 1) * 2);
    int tot = 0;

    set<int> is_in_hole;


    for (int i = 0; i < len_line; i++)
    {
        int u = hole_line[i].first, v = hole_line[i].second;
        nxt[++tot] = head[u];
        head[u] = tot;
        to[tot] = v;
        nxt[++tot] = head[v];
        head[v] = tot;
        to[tot] = u;

        if (is_in_hole.find(u) == is_in_hole.end())
            is_in_hole.insert(u);
        if (is_in_hole.find(v) == is_in_hole.end())
            is_in_hole.insert(v);

    }

    int* is_vis = (int*)malloc(sizeof(int) * l);
    memset(is_vis, 0, sizeof(int) * l);

    vector<int> hole_s;
    int* is_hole = (int*)malloc(sizeof(int) * l);
    memset(is_hole, 0, sizeof(int) * l);

    for (int i = 0; i < l / 3; i++)
    {
        if (is_in_hole.find(i) != is_in_hole.end())
        {
            point_color.push_back(1);
        }
        else
        {
            point_color.push_back(0);
        }
    }
    free(head); free(to); free(nxt);

    s.clear();
    return ans;

}

bool cmp_nnangle(const nnangle& a, const nnangle& b)
{
    if (a.visible == b.visible)
        return a.angle < b.angle;
    return a.visible;
}
vector<point> get_plain_base(vector<point>& points, vector<point>& fa)
{
    int len = points.size();
    point near = points[1];
    list<triangle> ans;
    /*
    if (len <= 2)
        return ans;
        */
    list<point> othpoint;
    list<point> hullpoint;
    point cur_fa = fa[points[0].pos];
    cur_fa.norm();

    float mam = 0;
    for (int i = 0; i < 3; i++)
        mam = max(mam, fabs((near - points[0]).indices[i]));

    int mq_size = min(20, (int)points.size() - 1);
    priority_queue<pair<float, point>, vector<pair<float, point> >, pqcmp >q;
    for (int i = 1; i < len; i++)
    {
        float t = l2_distance(points[i], points[0]);
        if ((int)q.size() >= mq_size)
        {
            float ma = q.top().first;
            if (t < ma)
            {
                q.pop();
                q.push(make_pair(t, points[i]));
            }
        }
        else
            q.push(make_pair(t, points[i]));
    }
    point mq_point[20];
    int mq_tot = 0;
    while (!q.empty())
    {
        mq_point[mq_tot++] = q.top().second; q.pop();
    }
    point v1 = (near - points[0]) - fa[near.pos] * dot_mul(cur_fa, (near - points[0]));
    point vi[20];
    float kpvi[20];
    pq_node angles[20];
    int mq_pos = -1;
    point mq_ver;
    for (int i = 0; i < mq_size; i++)
    {
        point tmp = (mq_point[i] - points[0]);
        vi[i] = tmp - fa[mq_point[i].pos] * dot_mul(cur_fa, tmp);
        kpvi[i] = 2 * dot_mul(cur_fa - fa[mq_point[i].pos], tmp) / (tmp.length() * tmp.length());
        angles[i].dis = acos(dot_mul(vi[i], v1) / (vi[i].length() * v1.length()));
        angles[i].pos = i;
        if (mq_pos == -1 && angles[i].dis < pi / 2)
        {
            mq_pos = i;
            float check1 = v1.length() * v1.length();
            float check2 = dot_mul(vi[i], v1) / check1;
            mq_ver = vi[i] - v1 * check2;
            //mq_ver = vi[i] - v1 * dot_mul(vi[i], v1) / (v1.length() * v1.length());
            mq_ver.norm();
        }
    }
    for (int i = 0; i < mq_size; i++)
    {
        if (dot_mul(vi[i], mq_ver) < 0)
            angles[i].dis = 2 * (float)pi - angles[i].dis;
    }
    sort(angles, angles + mq_size);
    float mq[3][3] = { 0 };
    for (int i = 0; i < mq_size; i++)
    {
        float tmp = angles[(i + 1) % 20].dis - angles[(i + 20) % 20].dis;
        if (tmp < 0)
            tmp += 2 * float(pi);
        tmp = tmp / 4 / (float)pi * kpvi[angles[i].pos];
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
            {
                mq[j][k] += vi[angles[i].pos].indices[j] * vi[angles[i].pos].indices[k];
            }
    }
    Eigen::Matrix3f matrix;
    matrix << mq[0][0], mq[0][1], mq[0][2], mq[1][0], mq[1][1], mq[1][2], mq[2][0], mq[2][1], mq[2][2];
    Eigen::EigenSolver<Eigen::Matrix3f> es(matrix);
    Eigen::Matrix3f D = es.pseudoEigenvalueMatrix();
    Eigen::Matrix3f V = es.pseudoEigenvectors();
    vector<point> pc(3);
    float kc[3];
    int mi;
    float mim = 1e7;
    for (int i = 0; i < 3; i++)
    {
        if (D(i, i) < mim)
        {
            mim = D(i, i);
            mi = i;
        }
    }

    for (int tot = 0, i = 0; i < 3; i++)
    {
        if (i != mi)
        {
            kc[tot] = D(i, i);
            for (int j = 0; j < 3; j++)
            {
                pc[tot].indices[j] = V(j, i);
            }
            tot++;
        }
    }
    float ktmp = kc[0]; kc[0] = 3 * kc[1] - kc[0]; kc[1] = 3 * ktmp - kc[1];
    point ptmp = pc[0]; pc[0] = pc[1] * 3 - pc[0]; pc[1] = ptmp * 3 - pc[1];
    if (kc[0] > kc[1])
    {
        ktmp = kc[0]; kc[0] = kc[1]; kc[1] = ktmp;
        ptmp = pc[0]; pc[0] = pc[1]; pc[1] = ptmp;
    }
    pc[0].norm(); pc[1].norm();
    return pc;
}
bool is_visible(point check_point, point line_first, point line_second, point origin)
{
    float a0 = line_first.indices[1] - line_second.indices[1];
    float b0 = line_second.indices[0] - line_first.indices[0];
    float c0 = line_first.indices[0] * line_second.indices[1] - line_first.indices[1] * line_second.indices[0];
    float a1 = -check_point.indices[1];
    float b1 = check_point.indices[0];
    float c1 = 0;
    if (origin.indices[0] != 0 && origin.indices[1] != 0)
    {
        a1 += origin.indices[1];
        b1 -= origin.indices[0];
        c1 = origin.indices[0] * check_point.indices[1] - origin.indices[1] * check_point.indices[0];
    }
    float div = a0 * b1 - a1 * b0;
    float x = (b0 * c1 - b1 * c0) / div;
    float y = (a1 * c0 - a0 * c1) / div;
    bool intersection_outside_twopoint;
    if (origin.indices[0] == 0 && origin.indices[1] == 0)
    {
        if (check_point.indices[0] > 0)
        {
            intersection_outside_twopoint = (x <= 0) || (x >= check_point.indices[0]);
        }
        else if (check_point.indices[0] < 0)
        {
            intersection_outside_twopoint = (x >= 0) || (x <= check_point.indices[0]);
        }
        else if (check_point.indices[1] > 0)
        {
            intersection_outside_twopoint = (y <= 0) || (y >= check_point.indices[1]);
        }
        else
            intersection_outside_twopoint = true;
    }
    else
    {
        if (check_point.indices[0] > origin.indices[0])
        {
            intersection_outside_twopoint = (x <= origin.indices[0]) || (x >= check_point.indices[0]);
        }
        else if (check_point.indices[0] < origin.indices[0])
        {
            intersection_outside_twopoint = (x >= origin.indices[0]) || (x <= check_point.indices[0]);
        }
        else if (check_point.indices[1] > origin.indices[1])
        {
            intersection_outside_twopoint = (y <= origin.indices[1]) || (y >= check_point.indices[1]);
        }
        else
            intersection_outside_twopoint = true;
    }
    if (intersection_outside_twopoint)
        return true;
    if (line_first.indices[0] > line_second.indices[0])
        return (x <= line_second.indices[0]) || (x >= line_first.indices[0]);
    if (line_first.indices[0] < line_second.indices[0])
        return (x >= line_second.indices[0]) || (x <= line_first.indices[0]);
    if (line_first.indices[1] > line_second.indices[1])
        return (y <= line_second.indices[1]) || (y >= line_first.indices[1]);
    if (line_first.indices[1] < line_second.indices[1])
        return (y >= line_second.indices[1]) || (y <= line_first.indices[1]);
    return false;
}

