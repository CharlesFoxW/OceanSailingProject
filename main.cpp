#include <iostream>
#include "oceanSailing.h"

// camera parameters
double Theta = PI / 6;
double Phi = PI / 6;
double R = 1.5;

// mouse control
int g_vMousePos[2];
int g_iLeftMouseButton,g_iRightMouseButton;

int windowWidth, windowHeight;

OceanSailing myOceanSailing;

void showBoundingBox();

void myinit()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(90.0,1.0,0.01,1000.0);

    // set background color
    glClearColor(0.8, 0.8, 0.8, 0.0);

    glCullFace(GL_BACK);
    glEnable(GL_CULL_FACE);

    glShadeModel(GL_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    myOceanSailing = OceanSailing();

    return;
}

void reshape(int w, int h)
{
    // Prevent a divide by zero, when h is zero.
    // You can't make a window of zero height.
    if(h == 0)
        h = 1;

    glViewport(0, 0, w, h);

    // Reset the coordinate system before modifying
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // Set the perspective

    double aspectRatio = 1.0 * w / h;
    gluPerspective(60.0f, aspectRatio, 0.01f, 1000.0f);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    windowWidth = w;
    windowHeight = h;

    glutPostRedisplay();
}

void display() {

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // camera parameters are Phi, Theta, R
    gluLookAt(R * cos(Phi) * cos(Theta), R * sin(Phi) * cos(Theta), R * sin(Theta) + 0.8,
              0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

    showBoundingBox();


    // draw particles:
    myOceanSailing.drawScene();

    glutSwapBuffers();
}

void doIdle() {
    // TODO: The animation computation:
    myOceanSailing.alterSceneByForwardEuler(0.01);

    glutPostRedisplay();
}

void mouseMotion (int x, int y)
{
    g_vMousePos[0] = x;
    g_vMousePos[1] = y;
}

void mouseButton(int button, int state, int x, int y)
{
    switch (button)
    {
        case GLUT_LEFT_BUTTON:
            g_iLeftMouseButton = (state==GLUT_DOWN);
            break;
        case GLUT_RIGHT_BUTTON:
            g_iRightMouseButton = (state==GLUT_DOWN);
            break;
        default:
            break;
    }

    g_vMousePos[0] = x;
    g_vMousePos[1] = y;
}

/* converts mouse drags into information about rotation/translation/scaling */
void mouseMotionDrag(int x, int y)
{
    int vMouseDelta[2] = {x-g_vMousePos[0], y-g_vMousePos[1]};

    if (g_iRightMouseButton) // handle camera rotations
    {
        Phi += vMouseDelta[0] * 0.01;
        Theta += vMouseDelta[1] * 0.01;

        if (Phi>2*PI)
            Phi -= 2*PI;

        if (Phi<0)
            Phi += 2*PI;

        if (Theta>PI / 2 - 0.01) // dont let the point enter the north pole
            Theta = PI / 2 - 0.01;

        if (Theta<- PI / 2 + 0.01)
            Theta = -PI / 2 + 0.01;

        g_vMousePos[0] = x;
        g_vMousePos[1] = y;
    }
}

// gets called whenever a key is pressed
void keyboardFunc (unsigned char key, int x, int y) {
    switch (key) {
        case 27:    //ESC key
            exit(0);
            break;

        default: break;
    }
}


int main(int argc, char** argv) {

    glutInit(&argc,argv);

    /* double buffered window, use depth testing, 640x480 */
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    windowWidth = 800;
    windowHeight = 800;
    glutInitWindowSize (windowWidth, windowHeight);
    glutInitWindowPosition (100,0);
    glutCreateWindow ("Sailing On The Sea");

    /* tells glut to use a particular display function to redraw */
    glutDisplayFunc(display);

    /* replace with any animate code */
    glutIdleFunc(doIdle);

    /* callback for mouse drags */
    glutMotionFunc(mouseMotionDrag);

    /* callback for window size changes */
    glutReshapeFunc(reshape);

    /* callback for mouse movement */
    glutPassiveMotionFunc(mouseMotion);

    /* callback for mouse button changes */
    glutMouseFunc(mouseButton);

    /* register for keyboard events */
    glutKeyboardFunc(keyboardFunc);

    /* do initialization */
    myinit();

    /* forever sink in the black hole */
    glutMainLoop();


    return 0;
}

void showBoundingBox()
{
    // draw the ground:
    glFrontFace(GL_CCW);
    glBegin(GL_TRIANGLES);
    glColor4f(0.3,0.3,0.3,1.0);
    glVertex3f((float) 0, (float) 0, (float) -0.2);
    glVertex3f((float) -0.5, (float) 0.5, (float) -0.2);
    glVertex3f((float) -0.5, (float) -0.5, (float) -0.2);

    glColor4f(0.2,0.2,0.2,1.0);
    glVertex3f((float) 0, (float) 0, (float) -0.2);
    glVertex3f((float) -0.5, (float) -0.5, (float) -0.2);
    glVertex3f((float) 0.5, (float) -0.5, (float) -0.2);
    glEnd();

    glFrontFace(GL_CW);
    glBegin(GL_TRIANGLES);
    glColor4f(0.3,0.3,0.3,1.0);

    glVertex3f((float) 0, (float) 0, (float) -0.2);
    glVertex3f((float) 0.5, (float) 0.5, (float) -0.2);
    glVertex3f((float) 0.5, (float) -0.5, (float) -0.2);

    glColor4f(0.2,0.2,0.2,1.0);
    glVertex3f((float) 0, (float) 0, (float) -0.2);
    glVertex3f((float) -0.5, (float) 0.5, (float) -0.2);
    glVertex3f((float) 0.5, (float) 0.5, (float) -0.2);
    glEnd();

    glFrontFace(GL_CW);
    glBegin(GL_TRIANGLES);
    glColor4f(0.3,0.3,0.3,1.0);
    glVertex3f((float) 0, (float) 0, (float) -0.2);
    glVertex3f((float) -0.5, (float) 0.5, (float) -0.2);
    glVertex3f((float) -0.5, (float) -0.5, (float) -0.2);

    glColor4f(0.2,0.2,0.2,1.0);
    glVertex3f((float) 0, (float) 0, (float) -0.2);
    glVertex3f((float) -0.5, (float) -0.5, (float) -0.2);
    glVertex3f((float) 0.5, (float) -0.5, (float) -0.2);
    glEnd();

    glFrontFace(GL_CCW);
    glBegin(GL_TRIANGLES);
    glColor4f(0.3,0.3,0.3,1.0);

    glVertex3f((float) 0, (float) 0, (float) -0.2);
    glVertex3f((float) 0.5, (float) 0.5, (float) -0.2);
    glVertex3f((float) 0.5, (float) -0.5, (float) -0.2);

    glColor4f(0.2,0.2,0.2,1.0);
    glVertex3f((float) 0, (float) 0, (float) -0.2);
    glVertex3f((float) -0.5, (float) 0.5, (float) -0.2);
    glVertex3f((float) 0.5, (float) 0.5, (float) -0.2);
    glEnd();

    glColor4f(0.4,0.4,0.4,1.0);
    glBegin(GL_LINES);
    glVertex3f(0.5, 0.5, -0.2);
    glVertex3f(0.5, 0.5, 0.8);

    glVertex3f(-0.5, 0.5, -0.2);
    glVertex3f(-0.5, 0.5, 0.8);

    glVertex3f(0.5, -0.5, -0.2);
    glVertex3f(0.5, -0.5, 0.8);

    glVertex3f(-0.5, -0.5, -0.2);
    glVertex3f(-0.5, -0.5, 0.8);
    glEnd();
}