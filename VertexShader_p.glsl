#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 point_color;


out vec3 mColor;

uniform mat4 view;

uniform mat4 projection;

void main()
{

    
    gl_Position = projection * view * vec4(aPos, 1.0);

	gl_PointSize = 5.0f;
    mColor = point_color;
}