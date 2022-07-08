#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include "shader.h"
#include "Camera.h"
#include "util.h"
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
#include<omp.h>
#include<io.h>
#include<time.h>
#include<thread>
#include<mutex>

//#include "data.h"
//#include "onlinetree.h"
//#include "onlinerf.h"

// #include"greedy_triangle.h"

#define pi 3.1415926
// #define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_DEBUG 
using namespace std;

unsigned int SCR_WIDTH = 800;
unsigned int SCR_HEIGHT = 600;
int draw_pos = 0;
int nums;
int num_camera = 0;
int draw_camera = 0;
bool is_pause = true;
glm::vec3 lightPos(1.2f, 1.0f, 2.0f);
glm::mat4 model = glm::mat4(1.0f);
double rad = 0, theta = 0, phi = 0,yangle=0,xangle=0;
Camera camera(glm::vec3(0.0f,0.0f,20.0f));//20
Camera image_camera(glm::vec3(0.0f, 0.0f, 2.0f));
bool firstMouse = true,is_left=true;
double lastX = SCR_WIDTH / 2.0f;
double lastY = SCR_HEIGHT / 2.0f;
float deltaTime = 0.0f;	// time between current frame and last frame
float lastFrame = 0.0f;
float cur_frame = 0.0f;
int pmode = 0;
int pos_mode = 0;
int pos_side = 0;
float cur_camera_center_x = 0,cur_camera_center_y = 0;
float cur_camera_center_last_x = SCR_WIDTH / 2.0f, cur_camera_center_last_y = SCR_HEIGHT / 2.0f;
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    
    if (firstMouse)
    {
        return;
    }

    if (is_left)
    {
        double dx = -(xpos - lastX) / SCR_WIDTH * pi / 2, dy = -(ypos - lastY) / SCR_HEIGHT * pi / 2;
        yangle += dx;
        xangle += dy;
        yangle -= floor(yangle / 2 / pi) * 2 * pi;
        xangle -= floor(xangle / 2 / pi) * 2 * pi;
        //cout << xangle << ' ' << yangle << endl;
        lastX = xpos;
        lastY = ypos;
    }
    else
    {
        double dx = (xpos - lastX) / SCR_WIDTH * pi / 2, dy = -(ypos - lastY) / SCR_HEIGHT * pi / 2;
        phi += dy;
        theta += dx;
        double x = rad * sin(phi) * cos(theta), y = rad * cos(phi), z = rad * sin(phi) * sin(theta);
        //camera.Position = glm::vec3(x, y, z) + camera.center;
        //cout << x << ' ' << y << ' ' << z << endl;
        lightPos = glm::vec3(x, y, z) + camera.center;
        lastX = xpos;
        lastY = ypos;
        if(is_pause)
            draw_pos = int((dx * 2 / pi) * nums/2 + draw_pos + nums) % nums;
        
    }
}

void scroll_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (ypos > 0)
    {
        if(pos_mode == 0)
            camera.ProcessKeyboard(FORWARD, deltaTime);
        else
        {
            double nx, ny;
            glfwGetCursorPos(window, &nx, &ny);
            cur_camera_center_x = nx; cur_camera_center_y = ny;
            if (nx < SCR_WIDTH / 4)
                pos_side = 0;
            else
                pos_side = 1;
            image_camera.ProcessKeyboard(FORWARD, deltaTime);
        }
    }
    else
    {
        if(pos_mode == 0)
            camera.ProcessKeyboard(BACKWARD, deltaTime);
        else
        {
            double nx, ny;
            glfwGetCursorPos(window, &nx, &ny);
            cur_camera_center_x = nx; cur_camera_center_y = ny;
            if (nx < SCR_WIDTH / 4)
                pos_side = 0;
            else
                pos_side = 1;
            image_camera.ProcessKeyboard(BACKWARD, deltaTime);
        }
    }
        
}
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
    {
        if (firstMouse)
        {
            glfwGetCursorPos(window, &lastX, &lastY);
            firstMouse = false;
            is_left = false;
        }
    }
    else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE)
    {
        firstMouse = true;
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        if (firstMouse)
        {
            glfwGetCursorPos(window, &lastX, &lastY);
            firstMouse = false;
            is_left = true;
        }
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    {
        firstMouse = true;
    }
    else if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_PRESS)
    {
        //pmode = (pmode + 1) % 3;
        pos_mode = pos_mode ^ 1;
        //is_pause = !is_pause;
    }
}
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    
    if (is_pause && (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS || glfwGetKey(window,GLFW_KEY_LEFT) == GLFW_PRESS))
    {
        //camera.ProcessKeyboard(LEFT, deltaTime);
        cur_frame += deltaTime;
        if (cur_frame > 0.1)
        {
            draw_pos = (draw_pos+nums-1)% nums;
            cur_frame = 0;
        }
        cout << "cur frame:" << draw_pos << endl;
    }
    if (is_pause && (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS))
    {
        //camera.ProcessKeyboard(RIGHT, deltaTime);
        cur_frame += deltaTime;
        if (cur_frame > 0.1)
        {
            draw_pos = (draw_pos+1)%nums;
            cur_frame = 0;
        }
        cout << "cur frame:" << draw_pos << endl;
    }
    

}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)

{
    SCR_WIDTH = width;
    SCR_HEIGHT = height;
    //glViewport(0, 0, width, height);

}

void readFromFile(const char* filepath,vector<float>& vertices)
{
    //cout << filepath;
    ifstream fin(filepath);
    //cout << fin.is_open();
    string s;
    int num = 0,i = -1;
    while (getline(fin, s))
    {
        vector<string> st;
        int pos = 0;
        while ((pos = s.find(' ')) != s.npos)
        {

            st.push_back(s.substr(0, pos ));
            s = s.substr(pos+1,s.length()-pos-1);

        }
        st.push_back(s);
        if (st[0][0] == 'e' && st[0][1] == 'l' && st[1][0] == 'v' && st[1][1] == 'e')
        {
            num = stoi(st[2]);
            continue;
        }
        if (st[0][0] == 'e' && st[0][1] == 'n')
        {
            i = 0;
            continue;
        }
        if (i >= 0 && i < num)
        {
            vertices.push_back((float)atof(st[0].c_str()));
            vertices.push_back((float)atof(st[1].c_str()));
            vertices.push_back((float)atof(st[2].c_str()));
            i++;
        } 
    }
    fin.close();

}

vector<float> chen(vector<float>& p, vector<float>& q)
{
    vector<float> ans;
    ans.push_back(p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3]);
    ans.push_back(p[0] * q[1] + p[1] * q[0] - p[2] * q[3] + p[3] * q[2]);
    ans.push_back(p[0] * q[2] + p[1] * q[3] + p[2] * q[0] - p[3] * q[1]);
    ans.push_back(p[0] * q[3] - p[1] * q[2] + p[2] * q[1] + p[3] * q[0]);
    return ans;
}
void write_ply(vector<float>& vertice, vector<face> indice,string filename)
{
    ofstream filen(filename.c_str());
    int l1 = vertice.size(), l2 = indice.size();
    filen << "ply"<<endl<<"format ascii 1.0"<<endl;
    filen << "element vertex " << l1 / 3 << endl;
    filen << "property float x" << endl;
    filen << "property float y" << endl;
    filen << "property float z" << endl;
    filen << "element face " << l2 << endl;
    filen << "property list uchar int vertex_indices" << endl;
    filen << "end_header" << endl;
    for (int i = 0; i < l1 / 3; i++)
        filen << vertice[i * 3] << " " << vertice[i * 3 + 1] << " " << vertice[i * 3 + 2] << endl;
    for (int i = 0; i < l2; i++)
        filen << "3 " << indice[i].indices[0] << " " << indice[i].indices[1] << " " << indice[i].indices[2] << endl;
    filen.close();
}


struct label
{
    float dimension[3];
    float location[3];
    float rotation_y;
};

struct cam_to_vel
{
    Eigen::Matrix4f R;
};

struct cell
{
    point points[8];
    bool label;
    cell()
    {
        label = false;
    }
    cell operator = (const cell& a)
    {
        for (int i = 0; i < 8; i++)
            this->points[i] = a.points[i];
        this->label = a.label;
        return *this;
    }
};

class proj_camera
{
public:
    point pos;
    point proj;
    point up;
    Eigen::Matrix<float,3,4> inner_matrix;
    Eigen::Matrix4f outer_matrix;
    int resolution = 10;
    proj_camera(point pos, point proj, point up, cell& cur_cell)
    {
        point cell_center;
        point min_cell;
        for (int i = 0; i < 8; i++)
            cell_center = cell_center + cur_cell.points[i];
        cell_center = cell_center / 8;
        point x_base = proj - pos; x_base.norm();
        point z_base = up - pos; z_base.norm();
        point y_base = x_mul3d(x_base, z_base); y_base.norm();
        cell lift_cell;
        for (int i = 0; i < 8; i++)
            lift_cell.points[i] = cur_cell.points[i] - pos;
        for (int i = 0; i < 8; i++)
        {
            point tmp_point;
            tmp_point.indices[0] = dot_mul(x_base, lift_cell.points[i]);
            tmp_point.indices[1] = dot_mul(y_base, lift_cell.points[i]);
            tmp_point.indices[2] = dot_mul(z_base, lift_cell.points[i]);
            lift_cell.points[i] = tmp_point;
            
        }
        Eigen::Matrix4f cur_cell_matrix; int ind[] = { 0,3,5,6 };
        for(int i = 0;i < 4;i++)
            for (int j = 0; j < 4; j++)
            {
                if (j == 3)
                    cur_cell_matrix(j, i) = 1;
                else 
                    cur_cell_matrix(j, i) = cur_cell.points[ind[i]].indices[j];
            }
        Eigen::Matrix4f result_outer_matrix;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                if (j == 3)
                    result_outer_matrix(j, i) = 1;
                else
                    result_outer_matrix(j, i) = lift_cell.points[ind[i]].indices[j];
            }
        outer_matrix = cur_cell_matrix.colPivHouseholderQr().solve(result_outer_matrix);
        cell_center = cell_center - pos;
        point tmp_point;
        for (int i = 0; i < 3; i++)
        {
            tmp_point.indices[0] = dot_mul(x_base, cell_center);
            tmp_point.indices[1] = dot_mul(y_base, cell_center);
            tmp_point.indices[2] = dot_mul(z_base, cell_center);
            cell_center = tmp_point;
        }
        for (int i = 0; i < 8; i++)
        {
            if (i == 0)
                min_cell = lift_cell.points[i];
            else
            {
                for (int j = 0; j < 3; j++)
                    min_cell.indices[j] = min(min_cell.indices[j], lift_cell.points[i].indices[j]);
            }
        }
        float u1 = ((float)resolution - 1.0f) / 2;
        float f_x = u1 / (cell_center.indices[0]/cell_center.indices[2] - min_cell.indices[0]/min_cell.indices[2]);
        float u_0 = -f_x * (min_cell.indices[0] / min_cell.indices[2]);
        float f_y = u1 / (cell_center.indices[1] / cell_center.indices[2] - min_cell.indices[1] / min_cell.indices[2]);
        float v_0 = -f_y * (min_cell.indices[1] / min_cell.indices[2]);
        inner_matrix(0, 0) = f_x; inner_matrix(0, 2) = u_0;
        inner_matrix(1, 1) = f_y; inner_matrix(1, 2) = v_0;
        inner_matrix(2, 2) = 1;
    }

};



struct camera_ksai
{
    Eigen::Vector3f trans;
    Eigen::Matrix3f rotation;
    int reference;
    int pos;
};

struct camera_display
{
    point points[5];
};
void read_slam_camera(const char* correct, const char* mine, vector<float*>& vertices, vector<unsigned int*>& indices,vector<int>& lens,vector<point>& centers,int& nums)
{
    point tmp_center = point();
    ifstream fin(correct);
    string s;
    vector<camera_ksai> correct_camera;
    vector<camera_ksai> mine_camera;
    vector<point> tmp_point;
    int num = 0, i = -1;
    while (getline(fin, s))
    {
        num++;
        vector<string> st;
        int pos = 0;

        while ((pos = s.find(' ')) != s.npos)
        {

            st.push_back(s.substr(0, pos));
            s = s.substr(pos + 1, s.length() - pos - 1);
        }
        if (s.length() != 0)
            st.push_back(s);
        camera_ksai tmp_camera;
        tmp_camera.reference = atoi(st[0].c_str());
        tmp_camera.pos = atoi(st[1].c_str());
        for (int j = 0; j < 3; j++)
            tmp_camera.trans[j] = atof(st[j + 2].c_str());
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                tmp_camera.rotation(j, k) = atof(st[5 + j * 3 + k].c_str());
        correct_camera.push_back(tmp_camera);
    }
    fin.close();
    ifstream fin2(mine);
    while (getline(fin2, s))
    {
        vector<string> st;
        int pos = 0;
        while ((pos = s.find(' ')) != s.npos)
        {

            st.push_back(s.substr(0, pos));
            s = s.substr(pos + 1, s.length() - pos - 1);
        }
        if (s.length() != 0)
            st.push_back(s);
        camera_ksai tmp_camera;
        tmp_camera.reference = atoi(st[0].c_str());
        tmp_camera.pos = atoi(st[1].c_str());
        for (int j = 0; j < 3; j++)
            tmp_camera.trans[j] = atof(st[j + 2].c_str());
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                tmp_camera.rotation(j, k) = atof(st[5 + j * 3 + k].c_str());
        mine_camera.push_back(tmp_camera);
    }
    fin2.close();
    tmp_point.push_back(point());
    nums = num;
    for (int j = 0; j < num; j++)
    {
        camera_ksai tmp_camera = correct_camera[j];
        Eigen::Vector3f ref = tmp_point[tmp_camera.reference].to_vector();
        Eigen::Vector3f tmp_ce = tmp_camera.rotation * ref + tmp_camera.trans;
        tmp_point.push_back(point(tmp_ce));
        tmp_center = tmp_center+ point(tmp_ce);
    }
    tmp_center = tmp_center / ((float)num+0.5);
    centers.push_back(tmp_center);
    float length = tmp_center.length() / ((float)num+0.5) *10;
    int flag = 0; float cmp_flag = 0;
    for (int i = 0; i < 3; i++)
    {
        if (fabs(tmp_center.indices[i]) > cmp_flag)
        {
            int f = 1;
            if (tmp_center.indices[i] < 0) 
                f = -1;
            flag = f * (i+1);
        }
    }
    camera_display initial;
    if (fabs(flag) == 0)
    {
        initial.points[4] = point(-1,length * flag, 0, 0);
        initial.points[0] = point(-1, 0, length, length);
        initial.points[1] = point(-1, 0, length,-length);
        initial.points[2] = point(-1, 0, -length , -length );
        initial.points[3] = point(-1, 0, -length , length);
    }
    else if (fabs(flag) == 1)
    {
        initial.points[4] = point(-1, 0, length * flag, 0);
        initial.points[0] = point(-1, length , 0, length );
        initial.points[1] = point(-1, length , 0, -length );
        initial.points[2] = point(-1, -length , 0, -length );
        initial.points[3] = point(-1, -length , 0, length );
    }
    else
    {
        initial.points[4] = point(-1, 0,0, length * flag);
        initial.points[0] = point(-1, length ,  length ,0);
        initial.points[1] = point(-1, length ,  -length ,0);
        initial.points[2] = point(-1, -length ,  -length ,0);
        initial.points[3] = point(-1, -length ,  length,0);
    }
    vector<camera_display> correct_displays;
    correct_displays.push_back(initial);
    for (int i = 0; i < num; i++)
    {
        camera_ksai tmp_camera = correct_camera[i];
        camera_display tmp_display;
        for (int j = 0; j < 5; j++)
        {
            Eigen::Vector3f ref = correct_displays[tmp_camera.reference].points[j].to_vector();
            Eigen::Vector3f tmp_ce = tmp_camera.rotation * ref + tmp_camera.trans;
            tmp_display.points[j] = point(tmp_ce);
        }
        correct_displays.push_back(tmp_display);
    }
    vector<camera_display> mine_displays;
    mine_displays.push_back(initial);
    for (int i = 0; i < num; i++)
    {
        camera_ksai tmp_camera = mine_camera[i];
        camera_display tmp_display;
        for (int j = 0; j < 5; j++)
        {
            Eigen::Vector3f ref = mine_displays[tmp_camera.reference].points[j].to_vector();
            Eigen::Vector3f tmp_ce = tmp_camera.rotation * ref + tmp_camera.trans;
            tmp_display.points[j] = point(tmp_ce);
        }
        mine_displays.push_back(tmp_display);
    }
    for (int i = 0; i < num; i++)
    {
        float* tmp_vertice = new float[(i + 1) * 10 * 9];
        for (int j = 0; j <= i; j++)
        {
            for (int k = 0; k < 5; k++)
            {
                int pos = j * 10 + k;
                for (int l = 0; l < 3; l++)
                    tmp_vertice[pos * 9 + l] = correct_displays[j].points[k].indices[l]*10;
                for (int l = 0; l < 3; l++)
                    tmp_vertice[pos * 9 + l + 3] = 0;
                tmp_vertice[pos * 9 + 6] = 1.0f;
                tmp_vertice[pos * 9 + 7] = 0.0f;
                tmp_vertice[pos * 9 + 8] = 0.0f;
            }
            for (int k = 0; k < 5; k++)
            {
                int pos = j * 10 + k + 5;
                for (int l = 0; l < 3; l++)
                    tmp_vertice[pos * 9 + l] = mine_displays[j].points[k].indices[l]*10;
                for (int l = 0; l < 3; l++)
                    tmp_vertice[pos * 9 + l + 3] = 0;
                tmp_vertice[pos * 9 + 6] = 0.0f;
                tmp_vertice[pos * 9 + 7] = 0.0f;
                tmp_vertice[pos * 9 + 8] = 1.0f;
            }
        }
        vertices.push_back(tmp_vertice);
        lens.push_back((i + 1) * 10);

        unsigned int* tmp_indices = new unsigned int[(i + 1) * 24];
        for (int j = 0; j <= i; j++)
        {
            int initial_pos = j * 10;
            for (int k = 0; k < 4; k++)
            {
                tmp_indices[j * 24 + k * 3] = initial_pos + k%4;
                tmp_indices[j * 24 + k * 3 + 1] = initial_pos + (k + 1)%4;
                tmp_indices[j * 24 + k * 3 + 2] = initial_pos + 4;
            }
            initial_pos += 5;
            for (int k = 0; k < 4; k++)
            {
                tmp_indices[j * 24 + k * 3+12] = initial_pos + k % 4;
                tmp_indices[j * 24 + k * 3 +12+ 1] = initial_pos + (k + 1) % 4;
                tmp_indices[j * 24 + k * 3 +12+ 2] = initial_pos + 4;
            }

        }
        indices.push_back(tmp_indices);

    }
}

void read_two_map_pointcloud(const char* correct, const char* mine, vector<float*>& vertices, vector<unsigned int*>& indices, vector<int>& lens, vector<point>& centers, int& nums)
{
    point tmp_center = point();
    ifstream fin(correct);
    string s;
    vector<point> actual_point;
    vector<point> my_point;
    int num = 0, i = -1;
    while (getline(fin, s))
    {
        vector<string> st;
        int pos = 0;

        while ((pos = s.find(' ')) != s.npos)
        {

            st.push_back(s.substr(0, pos));
            s = s.substr(pos + 1, s.length() - pos - 1);
        }
        if (s.length() != 0)
            st.push_back(s);
        if (atof(st[2].c_str()) == 0)
            continue;
        point tmp_point = point(-1,atof(st[0].c_str()), -atof(st[1].c_str()), -atof(st[2].c_str()));
        tmp_center = tmp_center + tmp_point;
        actual_point.push_back(tmp_point);
        num++;
    }
    ifstream myfin(mine);
    while (getline(myfin, s))
    {
        vector<string> st;
        int pos = 0;

        while ((pos = s.find(' ')) != s.npos)
        {

            st.push_back(s.substr(0, pos));
            s = s.substr(pos + 1, s.length() - pos - 1);
        }
        if (s.length() != 0)
            st.push_back(s);
        if (atof(st[2].c_str()) == 0)
            continue;
        point tmp_point = point(-1, atof(st[0].c_str()), -atof(st[1].c_str()), -atof(st[2].c_str()));
        tmp_center = tmp_center + tmp_point;
        my_point.push_back(tmp_point);
    }
    int l1 = actual_point.size();
    tmp_center = tmp_center / (l1 * 2);
    centers.push_back(tmp_center);
    nums = 1;
    float* vertice = new float[l1 * 2 * 9];
    for (int i = 0; i < l1; i++)
    {
        int pos = i * 9;
        for (int j = 0; j < 3; j++)
            vertice[pos + j] = actual_point[i].indices[j];
        for (int j = 3; j < 6; j++)
            vertice[pos + j] = 0;
        vertice[pos + 6] = 1.0f;
        vertice[pos + 7] = 0.0f;
        vertice[pos + 8] = 0.0f;
    }
    for (int i = 0; i < l1; i++)
    {
        int pos = i * 9 + l1 * 9;
        for (int j = 0; j < 3; j++)
            vertice[pos + j] = actual_point[i].indices[j];
        for (int j = 3; j < 6; j++)
            vertice[pos + j] = 0;
        vertice[pos + 6] = 0.0f;
        vertice[pos + 7] = 0.0f;
        vertice[pos + 8] = 0.0f;
    }
    vertices.push_back(vertice);
    lens.push_back(l1 * 2);
    unsigned int* indice;
    indices.push_back(indice);

}


void read_map_pointcloud(vector<float*>& vertices, vector<unsigned int*>& indices, vector<int>& lens, vector<point>& centers, int& nums)
{
    nums = 31;
    for (int i = 1; i <= nums; i++)
    {
        string filename = "D:\\lsd_slam\\cmp\\visualize\\" + to_string(i) + "_" + "my_point_var.txt";
        ifstream fin(filename.c_str());
        point tmp_center = point();
        string s;
        vector<point> my_point;

        while (getline(fin, s))
        {
            vector<string> st;
            int pos = 0;

            while ((pos = s.find(' ')) != s.npos)
            {

                st.push_back(s.substr(0, pos));
                s = s.substr(pos + 1, s.length() - pos - 1);
            }
            if (s.length() != 0)
                st.push_back(s);
            if (atof(st[2].c_str()) == 0)
                continue;
            point tmp_point = point(-1, atof(st[0].c_str()), -atof(st[1].c_str()), -atof(st[2].c_str()));
            tmp_center = tmp_center + tmp_point;

            my_point.push_back(tmp_point);
        }

        int len = my_point.size();
        tmp_center = tmp_center / len;
        centers.push_back(tmp_center);
        float* vertice = new float[len * 9];
        for (int j = 0; j < len; j++)
        {
            int pos = j * 9;
            for (int k = 0; k < 3; k++)
                vertice[pos + k] = my_point[j].indices[k];
            for (int k = 3; k < 6; k++)
                vertice[pos + k] = 0;
            vertice[pos + 6] = 0.0f;
            vertice[pos + 7] = 0.0f;
            vertice[pos + 8] = 0.0f;
        }
        vertices.push_back(vertice);
        lens.push_back(len);
        unsigned int* indice;
        indices.push_back(indice);
    }
}

void read_opt_camera(vector<float*>& vertices, vector<unsigned int*>& indices, vector<int>& num_points, vector<int>& num_indices, vector<point>& centers, int& nums)
{
    struct edge {
        int points[2];
    };
    vector<map<int, int> > idx2camera;
    string edge_file = "D:\\lsd_slam\\cmp\\sim3\\my_constraint.txt";
    string camera_file = "D:\\lsd_slam\\cmp\\sim3\\ac_camera.txt";
    ifstream edge_fin(edge_file.c_str());
    ifstream camera_fin(camera_file.c_str());
    string s;
    int num = 0;
    vector<vector<camera_ksai> > global_cameras;
    vector<vector<edge> > global_edges;
    vector<int> global_num_cameras;
    vector<int> global_num_edges;
    while (getline(camera_fin, s))
    {
        vector<camera_ksai> tmp_cameras;
        map<int, int> tmp_idx2camera;
        int num_camera = 0;
        if (s[0] >= '0' && s[0] <= '9')
            num_camera = atoi(s.c_str());
        if (num_camera == 0)
            continue;

        global_num_cameras.push_back(num_camera);
        point tmp_center;
        for (int i = 0; i < num_camera; i++)
        {
            string str;
            getline(camera_fin, str);
            vector<string> st;
            int pos = 0;

            while ((pos = str.find(' ')) != str.npos)
            {

                st.push_back(str.substr(0, pos));
                str = str.substr(pos + 1, str.length() - pos - 1);
            }
            if (str.length() != 0)
                st.push_back(str);
            camera_ksai tmp_ksai;
            tmp_ksai.pos = atoi(st[0].c_str());
            tmp_ksai.trans << atof(st[1].c_str()), atof(st[2].c_str()), atof(st[3].c_str());
            tmp_center = tmp_center + tmp_ksai.trans;
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    tmp_ksai.rotation(j, k) = atof(st[j * 3 + k + 4].c_str());
            tmp_cameras.push_back(tmp_ksai);
            tmp_idx2camera[tmp_ksai.pos] = i;
        }

        tmp_center = tmp_center / (tmp_cameras.size());
        centers.push_back(tmp_center);

        global_cameras.push_back(tmp_cameras);
        idx2camera.push_back(tmp_idx2camera);
        num++;
        string str;
        vector<edge> tmp_edges;
        while (getline(edge_fin, str))
        {
            if (str[0] == 'e')
                break;
            edge tmp_edge;
            vector<string> st;
            int pos = 0;

            while ((pos = str.find(' ')) != str.npos)
            {

                st.push_back(str.substr(0, pos));
                str = str.substr(pos + 1, str.length() - pos - 1);
            }
            if (str.length() != 0)
                st.push_back(str);
            int u = atoi(st[0].c_str()), v = atoi(st[1].c_str());
            if (u > v)
            {
                int t = u; u = v; v = t;
            }
            tmp_edge.points[0] = u; tmp_edge.points[1] = v;
            tmp_edges.push_back(tmp_edge);
            getline(edge_fin, str);
        }
        if (!global_edges.empty())
        {
            vector<edge> pre_edges = global_edges[global_edges.size() - 1];
            for (auto edge : pre_edges)
            {
                tmp_edges.push_back(edge);
            }


        }
        global_edges.push_back(tmp_edges);
        global_num_edges.push_back(tmp_edges.size());
    }

    nums = num;
    Eigen::Matrix3f K;
    K.setZero();
    K(0, 0) = 254.32695; K(1, 1) = 375.934387; K(0, 2) = 266.881897; K(1, 2) = 231.099091; K(2, 2) = 1;
    Eigen::Vector3f cameras[5];
    float w = 640 / 2/100.0, h = 480 / 2/100.0;
    cameras[0].setZero();
    for (int i = 1; i <= 4; i++)
    {
        cameras[i][0] = (i < 2|| i > 3) ? -w : w;
        cameras[i][1] = (i < 3) ? -h : h;
        cameras[i][2] = 1;
        cameras[i] = K.inverse() * cameras[i];
    }
    for (int i = 0; i < num;i++)
    {
        int num_point = global_num_cameras[i];
        float* vertice = new float[num_point * 5*9];
        for (int j = 0; j < num_point; j++)
        {
            for (int k = 0; k < 5; k++)
            {
                Eigen::Vector3f tmp_point = global_cameras[i][j].rotation * cameras[k] + global_cameras[i][j].trans;
                for (int l = 0; l < 3; l++)
                    vertice[j * 5 * 9 + k * 9 + l] = tmp_point[l];
                for (int l = 3; l < 6; l++)
                    vertice[j * 5 * 9 + k * 9 + l] = 0.0f;
                if (k == 0)
                {
                    vertice[j * 5 * 9 + k * 9 + 6] = 1.0f;
                    vertice[j * 5 * 9 + k * 9 + 7] = 0.0f;
                    vertice[j * 5 * 9 + k * 9 + 8] = 0.0f;
                }
                else
                {
                    vertice[j * 5 * 9 + k * 9 + 6] = 0.0f;
                    vertice[j * 5 * 9 + k * 9 + 7] = 0.0f;
                    vertice[j * 5 * 9 + k * 9 + 8] = 1.0f;
                }
            }
        }
        int num_edge = global_num_edges[i];
        unsigned int* indice = new unsigned int[num_edge * 2 + num_point * 8*2];
        for (int j = 0; j < num_edge; j++)
        {
            edge tmp_edge = global_edges[i][j];
            indice[j * 2] = idx2camera[i][tmp_edge.points[0]]*5;
            indice[j * 2 + 1] = idx2camera[i][tmp_edge.points[1]]*5;
        }
        for (int j = 0; j < num_point; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                for (int l = 0; l < 2; l++)
                {
                    indice[num_edge * 2 + j * 8 * 2 + k * 2 + l] = j * 5+l*(k+1);
                }
            }
            for (int k = 4; k < 8; k++)
            {
                for (int l = 0; l < 2; l++)
                {
                    indice[num_edge * 2 + j * 8 * 2 + k * 2 + l] = j * 5 + ((k - 4 + l) % 4)+1;
                }
            }
        }

        num_points.push_back(num_point*5);
        num_indices.push_back((num_edge + num_point * 8));

        vertices.push_back(vertice);
        indices.push_back(indice);
    }
    
}

void read_split(string& str,vector<string>& st)
{
    int pos = 0;
    st.clear();

    while ((pos = str.find(' ')) != str.npos)
    {

        st.push_back(str.substr(0, pos));
        str = str.substr(pos + 1, str.length() - pos - 1);
    }
    if (str.length() != 0)
        st.push_back(str);
}

void read_global_camera(vector<float*>& vertices, vector<unsigned int*>& indices, vector<int>& num_points, vector<int>& num_indices, int& nums)
{
    string camera_files = "D:\\lsd_slam\\cmp\\visualize\\pose.txt";
    ifstream fin(camera_files.c_str());
    string s;
    nums = 31;
    Eigen::Matrix3f K;
    K.setZero();
    K(0, 0) = 254.32695; K(1, 1) = 375.934387; K(0, 2) = 266.881897; K(1, 2) = 231.099091; K(2, 2) = 1;
    Eigen::Vector3f cameras[5];
    float w = 640 / 2 , h = 480 / 2;
    cameras[0].setZero();
    for (int i = 1; i <= 4; i++)
    {
        cameras[i][0] = (i < 2 || i > 3) ? -w : w;
        cameras[i][1] = (i < 3) ? -h : h;
        cameras[i][2] = -1;
        cameras[i] = K.inverse() * cameras[i];
    }
    for (int i = 0; i < nums; i++)
    {
        string str;
        getline(fin, str);
        vector<string> st;
        read_split(str, st);
        int u, v;
        camera_ksai camera_uv[2];
        u = atoi(st[0].c_str());
        v = atoi(st[1].c_str());
        string str1;
        getline(fin, str1);
        read_split(str1, st);
        camera_uv[0].pos = u;
        for (int j = 0; j < 3; j++)
            camera_uv[0].trans[j] = atof(st[j].c_str());
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                camera_uv[0].rotation(j, k) = atof(st[k+3*j+3].c_str());
        camera_uv[1].pos = v;
        string str2;
        getline(fin, str2);
        read_split(str2, st);
        for (int j = 0; j < 3; j++)
            camera_uv[1].trans[j] = atof(st[j].c_str());
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                camera_uv[1].rotation(j, k) = atof(st[k + 3 * j + 3].c_str());
        float* vertice = new float[2 * 5 * 9];
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k < 5; k++)
            {
                Eigen::Vector3f tmp_point = camera_uv[j].rotation.inverse() * cameras[k] - camera_uv[j].trans;
                for (int l = 0; l < 3; l++)
                    vertice[j * 5 * 9 + k * 9 + l] = tmp_point[l];
                for (int l = 3; l < 6; l++)
                    vertice[j * 5 * 9 + k * 9 + l] = 0.0f;
                if (j == 0)
                {
                    vertice[j * 5 * 9 + k * 9 + 6] = 1.0f;
                    vertice[j * 5 * 9 + k * 9 + 7] = 0.0f;
                    vertice[j * 5 * 9 + k * 9 + 8] = 0.0f;
                }
                else
                {
                    if (i == 13)
                    {
                        vertice[j * 5 * 9 + k * 9 + 6] = 1.0f;
                        vertice[j * 5 * 9 + k * 9 + 7] = 0.0f;
                        vertice[j * 5 * 9 + k * 9 + 8] = 0.0f;
                        continue;
                    }
                    vertice[j * 5 * 9 + k * 9 + 6] = 0.0f;
                    vertice[j * 5 * 9 + k * 9 + 7] = 0.0f;
                    vertice[j * 5 * 9 + k * 9 + 8] = 1.0f;
                }
            }
        }
        unsigned int* indice = new unsigned int[2 * 8 * 2];
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                for (int l = 0; l < 2; l++)
                {
                    indice[j * 8 * 2 + k * 2 + l] = j * 5 + l * (k + 1);
                }
            }
            for (int k = 4; k < 8; k++)
            {
                for (int l = 0; l < 2; l++)
                {
                    indice[j * 8 * 2 + k * 2 + l] = j * 5 + ((k - 4 + l) % 4) + 1;
                }
            }
        }

        num_points.push_back(2 * 5);
        num_indices.push_back((2 * 8));

        vertices.push_back(vertice);
        indices.push_back(indice);


        
    }

}
int u=3, v=3;
mutex uv_mutex;
void get_uv()
{
    int tmp_u, tmp_v;
    while (scanf("%d%d", &tmp_u, &tmp_v))
    {
        uv_mutex.lock();
        u = tmp_u;
        v = tmp_v;
        uv_mutex.unlock();
    }
}
int main()
{
    thread t1(get_uv);
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "OpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    Shader myShader("C:\\Users\\g107904\\source\\repos\\Project1\\Project1\\VertexShader.glsl", "C:\\Users\\g107904\\source\\repos\\Project1\\Project1\\FragmentShader.glsl");

    Shader imageShader("C:\\Users\\g107904\\source\\repos\\Project1\\Project1\\VertexShader2.glsl", "C:\\Users\\g107904\\source\\repos\\Project1\\Project1\\FragmentShader2.glsl");

    Shader pointShader("C:\\Users\\g107904\\source\\repos\\Project1\\Project1\\VertexShader_p.glsl", "C:\\Users\\g107904\\source\\repos\\Project1\\Project1\\FragmentShader_p.glsl");

    Shader cameraShader("C:\\Users\\g107904\\source\\repos\\Project1\\Project1\\VertexShader.glsl", "C:\\Users\\g107904\\source\\repos\\Project1\\Project1\\FragmentShader.glsl");
    //glEnable(GL_DEPTH_TEST);

    vector<string> files;
    //read_file_oder("C:\\Users\\g107904\\lyft_kitti\\order.txt", files);


    vector<float*> vertice;
    vector<point> fa;

    vector<point> centers;
    vector<int> lens;
    //get_crop_vertice(files);
    //read_all(files, centers, vertice,lens);

    //read_sample_data(files, centers, vertice, lens);
    

    //readFromFile("C:\\Users\\g107904\\Downloads\\bunny.tar\\bunny\\data\\bun000.ply", vertice);
    //collect_all_point(vertice);
    
    //readFromFile("C:\\Users\\g107904\\Downloads\\bunny.tar\\bunny\\reconstruction\\bun_zipper.ply", vertice);

   

    //vector<face> ans = greedy_triangle(vertice,fa,point_color);
    //vector<face> ans = pcl_greedy_triangle(vertice, fa, point_color);
    //vector<face> ans = grid_projection(vertice, fa, point_color);


    //vector<face> ans;

    //sort(ans.begin(), ans.end(), cmpind);

    /*
    for (int i = 0,pos = ini_l/3; i < len_label; i++,pos += 8)
    {
        ans.push_back(face(pos, pos + 1, pos + 2));
        ans.push_back(face(pos + 1, pos + 2, pos + 3));
        ans.push_back(face(pos + 4, pos + 5, pos + 6));
        ans.push_back(face(pos + 5, pos + 6, pos + 7));
        ans.push_back(face(pos, pos + 4, pos + 1));
        ans.push_back(face(pos + 4, pos + 1, pos + 5));
        ans.push_back(face(pos + 2, pos + 6, pos + 3));
        ans.push_back(face(pos + 6, pos + 3, pos + 7));
        ans.push_back(face(pos, pos + 2, pos + 4));
        ans.push_back(face(pos + 2, pos + 4, pos + 6));
        ans.push_back(face(pos + 1, pos + 5, pos + 3));
        ans.push_back(face(pos + 5, pos + 3, pos + 7));

    }
    */

    /*
    unsigned int* indices = new unsigned int[ans.size() * 3];
    for (int i = 0; i < (int)ans.size(); i++)
        for (int k = 0; k < 3; k++)
            indices[i * 3 + k] = ans[i].indices[k];
    */

    vector<unsigned int*> indices;

    /*
    unsigned int indices[] = {  // note that we start from 0!
        0, 1, 3,  // first Triangle
        1, 2, 3   // second Triangle
    };
    */

    //write_ply(vertice, ans);


    //cout << "face number" << ans.size() << endl;

    //nums = files.size();

    //read_slam_camera("D:\\lsd_slam\\cmp\\actual.txt", "D:\\lsd_slam\\cmp\\my.txt", vertice, indices,lens,centers,nums);

    //read_two_map_pointcloud("D:\\lsd_slam\\cmp\\map\\pointcloud\\actual_pointcloud.txt", "D:\\lsd_slam\\cmp\\map\\pointcloud\\my_pointcloud.txt", vertice, indices, lens, centers, nums);

    read_map_pointcloud(vertice, indices, lens, centers, nums);

    vector<int> num_points;
    vector<int> num_indices;
    //read_opt_camera(vertice, indices, num_points,num_indices, centers, nums);

    unsigned int* VBO = new unsigned int[nums];
    unsigned int* VAO = new unsigned int[nums];
    unsigned int* EBO = new unsigned int[nums];
    glGenVertexArrays(nums, VAO);
    glGenBuffers(nums, VBO);
    glGenBuffers(nums, EBO);

    for (int i = 0; i < nums; i++)
    {
        glBindVertexArray(VAO[i]);
        glBindBuffer(GL_ARRAY_BUFFER, VBO[i]);

        glBufferData(GL_ARRAY_BUFFER, lens[i]*9*sizeof(float), vertice[i], GL_STATIC_DRAW);
        //glBufferData(GL_ARRAY_BUFFER, num_points[i] * 9 * sizeof(float), vertice[i], GL_STATIC_DRAW);

        //position
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        //normal
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(1); 

        //point_color
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
        glEnableVertexAttribArray(2);
    

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO[i]);

        //glBufferData(GL_ELEMENT_ARRAY_BUFFER, lens[i]/5*4 * 3 * sizeof(float), indices[i], GL_STATIC_DRAW);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 0, indices[i], GL_STATIC_DRAW);
        //glBufferData(GL_ELEMENT_ARRAY_BUFFER, num_indices[i]*2*sizeof(unsigned int), indices[i], GL_STATIC_DRAW);
    }
    
    
    
    
    //glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    //glBindBuffer(GL_ARRAY_BUFFER, 0);

    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);

    //glBufferData(GL_ELEMENT_ARRAY_BUFFER, ans.size()*3*sizeof(float), indices, GL_STATIC_DRAW);
    
    glEnable(GL_PROGRAM_POINT_SIZE);
    
    //image

    nums = 31;
    
    float distance = 0.35f;

    int width = 640, height = 480;

    float image_vertices[]= {
        // positions          // colors           // texture coords
         -distance+0.32f,  -0.24f, 0.0f,   1.0f, 0.0f, 0.0f,   1.0f, 1.0f, // top right
         -distance+0.32f, 0.24f, 0.0f,   0.0f, 1.0f, 0.0f,   1.0f, 0.0f, // bottom right
        -distance-0.32f, 0.24f, 0.0f,   0.0f, 0.0f, 1.0f,   0.0f, 0.0f, // bottom left
        -distance-0.32f, -0.24f, 0.0f,   1.0f, 1.0f, 0.0f,   0.0f, 1.0f  // top left 
    };
    unsigned int image_indices[] = {
    0, 1, 3, // first triangle
    1, 2, 3  // second triangle
    };
    unsigned int image_VBO, image_VAO, image_EBO;

    glGenVertexArrays(1, &image_VAO);
    glGenBuffers(1, &image_VBO);
    glGenBuffers(1, &image_EBO);

    glBindVertexArray(image_VAO);

    glBindBuffer(GL_ARRAY_BUFFER, image_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(image_vertices), image_vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, image_EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(image_indices), image_indices, GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    // texture coord attribute
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
    
    unsigned int* texture = new unsigned int[nums];
    glGenTextures(nums, texture);
    vector<unsigned char*> kf_data;
    for (int it = 0; it < nums; it++)
    {
        string filename = "D:\\lsd_slam\\cmp\\visualize\\" + to_string(it+1) + "_" + "image_kf.txt";
        FILE* fp = fopen(filename.c_str(), "r+");
        unsigned char* data = new unsigned char[width * height * 3];
        for (int j = 0; j < height; j++)
            for (int i = 0; i < width; i++)
            {
                float tmp_image;
                fscanf(fp, "%f", &tmp_image);
                int pos = i + j * width;
                data[pos * 3] = (unsigned char)tmp_image;
                data[pos * 3 + 2] = data[pos * 3 + 1] = data[pos * 3];
            }
        fclose(fp);
        kf_data.push_back(data);

    }
    for (int it = 0; it < nums; it++)
    {
        glBindTexture(GL_TEXTURE_2D, texture[it]); // all upcoming GL_TEXTURE_2D operations now have effect on this texture object
        // set the texture wrapping parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	// set texture wrapping to GL_REPEAT (default wrapping method)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        // set texture filtering parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

       
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, kf_data[it]);
        //glGenerateMipmap(GL_TEXTURE_2D);
        glBindVertexArray(0);
    }


    float image_ref_vertices[] = {
        // positions          // colors           // texture coords
         distance + 0.32f,  -0.24f, 0.0f,   1.0f, 0.0f, 0.0f,   1.0f, 1.0f, // top right
         distance + 0.32f, 0.24f, 0.0f,   0.0f, 1.0f, 0.0f,   1.0f, 0.0f, // bottom right
        distance - 0.32f, 0.24f, 0.0f,   0.0f, 0.0f, 1.0f,   0.0f, 0.0f, // bottom left
        distance - 0.32f, -0.24f, 0.0f,   1.0f, 1.0f, 0.0f,   0.0f, 1.0f  // top left 
    };

    unsigned int image_ref_VBO, image_ref_VAO, image_ref_EBO;

    glGenVertexArrays(1, &image_ref_VAO);
    glGenBuffers(1, &image_ref_VBO);
    glGenBuffers(1, &image_ref_EBO);

    glBindVertexArray(image_ref_VAO);

    glBindBuffer(GL_ARRAY_BUFFER, image_ref_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(image_ref_vertices), image_ref_vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, image_ref_EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(image_indices), image_indices, GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    // texture coord attribute
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    unsigned int* image_ref_texture = new unsigned int[nums];
    glGenTextures(nums, image_ref_texture);
    vector<unsigned char*> image_ref_data;
    for (int it = 0; it < nums; it++)
    {
        string filename = "D:\\lsd_slam\\cmp\\visualize\\" + to_string(it+1) + "_" + "image_ref.txt";
        FILE* fp = fopen(filename.c_str(), "r+");
        unsigned char* data = new unsigned char[width * height * 3];
        for (int j = 0; j < height; j++)
            for (int i = 0; i < width; i++)
            {
                float tmp_image;
                fscanf(fp, "%f", &tmp_image);
                int pos = i + j * width;
                data[pos * 3] = (unsigned char)tmp_image;
                data[pos * 3 + 2] = data[pos * 3 + 1] = data[pos * 3];
            }
        fclose(fp);
        image_ref_data.push_back(data);
    }
    for (int it = 0; it < nums; it++)
    {
        glBindTexture(GL_TEXTURE_2D, image_ref_texture[it]); // all upcoming GL_TEXTURE_2D operations now have effect on this texture object
        // set the texture wrapping parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	// set texture wrapping to GL_REPEAT (default wrapping method)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        // set texture filtering parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);


        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image_ref_data[it]);
    }
    float image_depth_vertices[] = {
        // positions          // colors           // texture coords
           0.32f,  -0.5f-0.24f, 0.0f,   1.0f, 0.0f, 0.0f,   1.0f, 1.0f, // top right
          0.32f, -0.5f+0.24f, 0.0f,   0.0f, 1.0f, 0.0f,   1.0f, 0.0f, // bottom right
         - 0.32f, -0.5f+0.24f, 0.0f,   0.0f, 0.0f, 1.0f,   0.0f, 0.0f, // bottom left
         - 0.32f, -0.5f-0.24f, 0.0f,   1.0f, 1.0f, 0.0f,   0.0f, 1.0f  // top left 
    };

    unsigned int image_depth_VBO, image_depth_VAO, image_depth_EBO;

    glGenVertexArrays(1, &image_depth_VAO);
    glGenBuffers(1, &image_depth_VBO);
    glGenBuffers(1, &image_depth_EBO);

    glBindVertexArray(image_depth_VAO);

    glBindBuffer(GL_ARRAY_BUFFER, image_depth_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(image_depth_vertices), image_depth_vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, image_depth_EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(image_indices), image_indices, GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    // texture coord attribute
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    unsigned int* image_depth_texture = new unsigned int[nums];
    glGenTextures(nums, image_depth_texture);
    vector<unsigned char*> image_depth_data;
    for (int it = 0; it < nums; it++)
    {
        string filename = "D:\\lsd_slam\\cmp\\visualize\\" + to_string(it+1) + "_" + "image_depth.txt";
        FILE* fp = fopen(filename.c_str(), "r+");
        unsigned char* data = new unsigned char[width * height * 3];
        vector<float> tmp_data;
        for(int j = 0;j < height;j++)
            for (int i = 0; i < width; i++)
            {
                float tmp_image;
                fscanf(fp, "%f", &tmp_image);
                tmp_data.push_back(tmp_image);
            }
        for (int j = 0; j < height; j++)
            for (int i = 0; i < width; i++)
            {

                int pos = i + j * width;
                float tmp_image = tmp_data[pos];
                if (tmp_image <= 0)
                {
                    data[pos * 3 + 2] = data[pos * 3 + 1] = data[pos * 3] = 255;
                }
                else
                {
                    float r = (0.0f - tmp_image) * 255 / 1.0f; if (r < 0) r = -r;
                    float g = (1.0f - tmp_image) * 255 / 1.0f; if (g < 0) g = -g;
                    float b = (2.0f - tmp_image) * 255 / 1.0f; if (b < 0) b = -b;
                    unsigned char rc = r < 0 ? 0 : (r > 255 ? 255 : r);
                    unsigned char gc = g < 0 ? 0 : (g > 255 ? 255 : g);
                    unsigned char bc = b < 0 ? 0 : (b > 255 ? 255 : b);
                    data[pos * 3] = 255 - rc;
                    data[pos * 3 + 1] = 255 - gc;
                    data[pos * 3 + 2] = 255 - bc;
                }
            }
        fclose(fp);
        image_depth_data.push_back(data);
    }
    for (int it = 0; it < nums; it++)
    {
        glBindTexture(GL_TEXTURE_2D, image_depth_texture[it]); // all upcoming GL_TEXTURE_2D operations now have effect on this texture object
        // set the texture wrapping parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	// set texture wrapping to GL_REPEAT (default wrapping method)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        // set texture filtering parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);


        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image_depth_data[it]);
    }
    

    //glBindVertexArray(0);

    int pre_pos = 0;

    int nrAttributes;
    glGetIntegerv(GL_MAX_VERTEX_ATTRIBS, &nrAttributes);
    std::cout << "Maximum nr of vertex attributes supported: " << nrAttributes << std::endl;

    
    unsigned int point_VBO, point_VAO;
    float* point_vertices = 0;

    glGenVertexArrays(1, &point_VAO);
    glGenBuffers(1, &point_VBO);



    unsigned int pc_VBO, pc_VAO;
    float* pc_vertices = 0;

    glGenVertexArrays(1, &pc_VAO);
    glGenBuffers(1, &pc_VBO);
    

    vector<unsigned int*> camera_indices;
    vector<float*> camera_vertices;
    read_global_camera(camera_vertices, camera_indices, num_points, num_indices, num_camera);

    unsigned int* camera_VBO = new unsigned int[num_camera];
    unsigned int* camera_VAO = new unsigned int[num_camera];
    unsigned int* camera_EBO = new unsigned int[num_camera];
    glGenVertexArrays(num_camera, camera_VAO);
    glGenBuffers(num_camera, camera_VBO);
    glGenBuffers(num_camera, camera_EBO);

    for (int i = 0; i < num_camera; i++)
    {
        glBindVertexArray(camera_VAO[i]);
        glBindBuffer(GL_ARRAY_BUFFER, camera_VBO[i]);

        glBufferData(GL_ARRAY_BUFFER, num_points[i] * 9 * sizeof(float), camera_vertices[i], GL_STATIC_DRAW);

        //position
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        //normal
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(1);

        //point_color
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
        glEnableVertexAttribArray(2);


        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO[i]);

        glBufferData(GL_ELEMENT_ARRAY_BUFFER, num_indices[i]*2*sizeof(unsigned int), camera_indices[i], GL_STATIC_DRAW);
    }



    long long prev_time = clock();

    point tmp_center;

    while (!glfwWindowShouldClose(window))
    {
        //point center = centers[draw_pos];

        uv_mutex.lock();
        int flag = 0;

        float currentFrame = (float)glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        processInput(window);

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

        //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glClear(GL_COLOR_BUFFER_BIT);        
        
        //image
        
        int cur_pos = u + v * 640;
        /*
        if (cur_pos != pre_pos)
        {
            pre_pos = cur_pos;
            flag = 1;
        }
        */
        image_camera.center = glm::vec3(tmp_center.indices[0], tmp_center.indices[1], tmp_center.indices[2]);
        glm::mat4 image_view = image_camera.GetViewMatrix();

        glm::mat4 image_projection = glm::perspective(glm::radians(image_camera.Zoom), 1.0f, 0.1f, 100.0f);

        //image_view = glm::mat4(1.0f);
        //image_projection = glm::mat4(1.0f);

        glViewport(0, 0, SCR_WIDTH / 2, SCR_HEIGHT);

        imageShader.use();

        imageShader.setMat4("view", image_view);
        imageShader.setMat4("projection", image_projection);
        
        glBindTexture(GL_TEXTURE_2D, texture[draw_pos]);
        glBindVertexArray(image_VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);


        glBindTexture(GL_TEXTURE_2D, image_depth_texture[draw_pos]);
        glBindVertexArray(image_depth_VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        
        glBindTexture(GL_TEXTURE_2D, image_ref_texture[draw_pos]);        
        glBindVertexArray(image_ref_VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        
        
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        


        pointShader.use();

        pointShader.setMat4("view", image_view);
        pointShader.setMat4("projection", image_projection);

        if (flag)
        {
            cout << "pixel value:" << (float)(*(kf_data[draw_pos] + (u + v * width) * 3)) << endl;
            point image_p = point(-1, (float)u / 640 * 0.64 + (-distance - 0.32f), -(float)v / 480 * 0.48 + 0.24f, 0.0f);
            point image_ref_p = point();
            point image_depth_p = point(-1, (float)u / 640 * 0.64 + (-0.32f), -(float)v / 480 * 0.48 + (-0.5f + 0.24f), 0.0f);

            string str = "D:\\lsd_slam\\cmp\\map\\"+to_string(draw_pos)+"_"+"map\\" + to_string(cur_pos) + ".txt";
            FILE* image_fp = fopen(str.c_str(),"r+");
            if (!image_fp)
                cout << "grad not satisfy!\n" << endl;
            else
            {
                float epl_x,epl_y;
                fscanf(image_fp, "%f", &epl_x);
                if (epl_x == 0)
                    cout << "tracking is not good at this point! or epl is not good at this point!\n" << endl;
                else
                {
                    fscanf(image_fp, "%f", &epl_y);
                    int is_good;
                    fscanf(image_fp, "%d", &is_good);
                    if (!is_good)
                    {
                        cout << "epl is not good!" << endl;
                    }
                    else
                    {
                        point kf_point[5], ref_point[5];
                        float match_x=0, match_y=0;
                        for (int j = 0; j < 5; j++)
                        {
                            fscanf(image_fp, "%f", &(kf_point[j].indices[0]));
                            fscanf(image_fp, "%f", &(kf_point[j].indices[1]));
                        }
                        point p_close, p_far;
                        fscanf(image_fp, "%f%f", &(p_close.indices[0]), &(p_close.indices[1]));
                        fscanf(image_fp, "%f%f", &(p_far.indices[0]), &(p_far.indices[1]));
                        for (int j = 0; j < 5; j++)
                        {
                            fscanf(image_fp, "%f", &(ref_point[j].indices[0]));
                            if (j == 0 && ref_point[j].indices[0] == 0)
                            {
                                is_good = 0;
                                break;
                            }
                            fscanf(image_fp, "%f", &(ref_point[j].indices[1]));
                        }
                        if (!is_good)
                        {
                            cout << "p_close or p_far not satisfy!" << endl;
                        }
                        for (int i = 0;is_good; i++)
                        {
                            int it;
                            fscanf(image_fp, "%d", &it);
                            if (it < i)
                            {
                                if (it == -1)
                                {
                                    cout << "find error too large!" << endl;
                                    break;
                                }
                                if (it == 0)
                                {
                                    cout << "did subpixel fail!" << endl;
                                    break;
                                }
                                else
                                {
                                    
                                    fscanf(image_fp, "%f%f", &match_x, &match_y);
                                    break;
                                }
                            }
                            fscanf(image_fp, "%f%f", &match_x, &match_y);
                            float np_x, np_y;
                            fscanf(image_fp, "%f%f", &np_x, &np_y);

                        }
                        cout << "matched:" << match_x << ' ' << match_y << endl;
                        
                        int match_px = (int)match_x, match_py = (int)match_y;
                        float match_a = match_x - match_px, match_b = match_y - match_py;
                        float tmp_pixel = 0;
                        tmp_pixel += match_a * match_b * (*(image_ref_data[draw_pos] + (match_px + match_py * width) * 3));
                        tmp_pixel += match_a * (1 - match_b) * (*(image_ref_data[draw_pos]+(match_px+(match_py+1)*width)*3));
                        tmp_pixel += (1 - match_a) * match_b * (*(image_ref_data[draw_pos] + (match_px+1 + match_py * width) * 3));
                        tmp_pixel += (1 - match_a) * (1 - match_b) * (*(image_ref_data[draw_pos] + (match_px + 1 + (match_py+1) * width) * 3));

                        cout << "match pixel value:" << tmp_pixel << endl;
                        image_ref_p = point(-1, (float)match_x / 640 * 0.64 + distance - 0.32f, -(float)match_y / 480 * 0.48 + 0.24f, 0.0f);
                    }
                }
                fclose(image_fp);
            }
            
            point_vertices = new float[3 * 6];
            for (int i = 0; i < 3; i++)
            {
                point_vertices[0 + i] = image_p.indices[i];
                point_vertices[6 + i] = image_ref_p.indices[i];
                point_vertices[12 + i] = image_depth_p.indices[i];
                point_vertices[i * 6 + 3] = 0.0f;
                point_vertices[i * 6 + 3 + 1] = 0.0f;
                point_vertices[i * 6 + 3 + 2] = 1.0f;
            }
            glBindVertexArray(point_VAO);

            glBindBuffer(GL_ARRAY_BUFFER, point_VBO);
            glBufferData(GL_ARRAY_BUFFER, 18*sizeof(float), point_vertices, GL_STATIC_DRAW);

            // position attribute
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(0);
            // color attribute
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
            glEnableVertexAttribArray(1);
        }
        //glEnable(GL_PROGRAM_POINT_SIZE);
        glBindVertexArray(point_VAO);
        glDrawArrays(GL_POINTS, 0, 3);


        
        
        
        glViewport(SCR_WIDTH / 2, 0, SCR_WIDTH / 2, SCR_HEIGHT);
        
        

        
        

        //center.indices[2] = -5.0f;
        
        long long cur_time = clock();
        
        draw_pos = (draw_pos + nums) % nums;
        
        if (!is_pause && (cur_time-prev_time) % 60 == 0)
        {
            draw_pos = (draw_pos + 1) % nums;prev_time = cur_time;
        }

        
        point center = centers[draw_pos];

        camera.center = glm::vec3(center.indices[0], center.indices[1], center.indices[2]);



        //glm::vec3 dr = camera.Position - camera.center;
        glm::vec3 dr = lightPos - camera.center;
        rad = sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
        phi = acos(dr.y / rad);
        theta = acos(dr.x / rad / sin(phi));

        model = glm::mat4(1.0f);
        glm::vec3 tmp = camera.Position - camera.center;

        

        myShader.use();
        glm::vec3 lightColor = glm::vec3(1.0f,1.0f,1.0f);
        glm::vec3 diffuseColor = lightColor * glm::vec3(0.5f); // decrease the influence
        glm::vec3 ambientColor = diffuseColor * glm::vec3(0.2f); // low influence
        myShader.setVec3("light.ambient", ambientColor);
        myShader.setVec3("light.diffuse", diffuseColor);
        myShader.setVec3("light.specular", 1.0f, 1.0f, 1.0f);

        // material properties

        //myShader.setVec3("material.ambient", 0.5f, 0.5f, 0.5f);
        //myShader.setVec3("material.diffuse", 0.5f, 0.5f, 0.5f);
        //myShader.setVec3("material.specular", 0.5f, 0.5f, 0.5f); // specular lighting doesn't have full effect on this object's material
        //myShader.setFloat("material.shininess", 32.0f);

        myShader.setVec3("material.ambient", 0.2125f,0.1275f,0.054f);
        myShader.setVec3("material.diffuse", 0.714f,0.4284f,0.18144f);
        myShader.setVec3("material.specular", 0.393548f,0.271906f,0.166721f); // specular lighting doesn't have full effect on this object's material
        myShader.setFloat("material.shininess", 25.6f);

        myShader.setVec3("light.position", lightPos);
        myShader.setVec3("viewPos", camera.Position);

        

        model = glm::translate(model, camera.center);
        //model = glm::rotate(model, float(xangle), glm::vec3(1.0f, 0.0f, 0.0f));
        //model = glm::rotate(model, float(yangle), glm::vec3(0.0f, 1.0f, 0.0f));
        glm::qua<float> q = glm::qua<float>((glm::vec3(float(xangle), float(yangle), 0.0f)));
        model =  model * glm::mat4_cast(q);
        model = glm::translate(model, -camera.center);

        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        myShader.setMat4("projection", projection);

        myShader.setMat4("model", model);

        
        glm::mat4 view = camera.GetViewMatrix();
        myShader.setMat4("view", view);

        glBindVertexArray(VAO[draw_pos]); 
        glDrawArrays(GL_POINTS, 0, lens[draw_pos]);
        //glDrawElements(GL_LINES, num_indices[draw_pos]*2, GL_UNSIGNED_INT, 0);
        //glDrawArrays(GL_POINTS, 0, num_points[draw_pos]);
        
        cameraShader.use();
        cameraShader.setVec3("light.ambient", ambientColor);
        cameraShader.setVec3("light.diffuse", diffuseColor);
        cameraShader.setVec3("light.specular", 1.0f, 1.0f, 1.0f);
        cameraShader.setVec3("material.ambient", 0.2125f, 0.1275f, 0.054f);
        cameraShader.setVec3("material.diffuse", 0.714f, 0.4284f, 0.18144f);
        cameraShader.setVec3("material.specular", 0.393548f, 0.271906f, 0.166721f); // specular lighting doesn't have full effect on this object's material
        cameraShader.setFloat("material.shininess", 25.6f);
        cameraShader.setVec3("light.position", lightPos);
        cameraShader.setVec3("viewPos", camera.Position);
        cameraShader.setMat4("projection", projection);
        cameraShader.setMat4("model", model);
        cameraShader.setMat4("view", view);

        glBindVertexArray(camera_VAO[draw_pos]);
        glDrawElements(GL_LINES, num_indices[draw_pos] * 2, GL_UNSIGNED_INT, 0);
        //glDrawArrays(GL_POINTS, 0, num_points[draw_pos]);





    /*
    if (flag)
    {
        float* cur_vertice = vertice[0];
        point tmp_point;
        pc_vertices = new float[9];
        for (int i = 0; i < 3; i++)
        {
            pc_vertices[i] = *(cur_vertice + (u * 480 + v) * 9 + i);
        }

        pc_vertices[3] = 0.0f;
        pc_vertices[4] = 0.0f;
        pc_vertices[5] = 0.0f;
        pc_vertices[6] = 0.0f;
        pc_vertices[7] = 0.0f;
        pc_vertices[8] = 1.0f;

        glBindVertexArray(pc_VAO);

        glBindBuffer(GL_ARRAY_BUFFER, pc_VBO);
        glBufferData(GL_ARRAY_BUFFER, 9 * sizeof(float), pc_vertices, GL_STATIC_DRAW);


        // position attribute
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(1);
        // color attribute
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
        glEnableVertexAttribArray(2);
    }
    glEnable(GL_PROGRAM_POINT_SIZE);
    glBindVertexArray(pc_VAO);
    glDrawArrays(GL_POINTS, 0, 1);
    */
        /*
        if (pmode == 0)
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
            glDrawArrays(GL_POINTS, 0, lens[draw_pos]);
        }
        else if (pmode == 1)
        {

            //glDepthMask(GL_FALSE);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDrawElements(GL_TRIANGLES, lens[draw_pos]/5*4 * 3, GL_UNSIGNED_INT, 0);

            //glDepthMask(GL_TRUE);
        }
        else
        {

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glDrawElements(GL_TRIANGLES, lens[draw_pos]/5*4 * 3, GL_UNSIGNED_INT, 0);
        }
        */
        /*
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        //glDrawArrays(GL_TRIANGLES, 0, l/3);
        glDrawElements(GL_TRIANGLES, ans.size()*3, GL_UNSIGNED_INT, 0);
        //glDrawArrays(GL_POINTS, 0, l);
        */
        flag = 0;
        uv_mutex.unlock();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();
    return 0;
}
