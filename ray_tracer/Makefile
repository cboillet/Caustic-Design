EXECUTABLE = cg1_ex5
CC = g++
CPPFLAGS = -O3 -Iglm/
DEPS = glm/glm.hpp  /glm/gtc/matrix_transform.hpp  /glm/gtc/matrix_inverse.hpp Window.hpp Context.hpp TriMesh.hpp GLSLShader.hpp Assets.hpp Scene.hpp View3D.hpp ViewRender.hpp Model.hpp Ray.hpp Image.hpp
OBJ = Window.o Context.o TriMesh.o GLSLShader.o Assets.o Scene.o View3D.o ViewRender.o Model.o Ray.o Image.o
LIBS = -lGLEW -lGLU -lGL -lglut

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CPPFLAGS)

$(EXECUTABLE): $(OBJ)
	$(CC) -o $@ $^ $(CPPFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o *~ $(EXECUTABLE)		
