#version 330 core
layout(location = 0) out vec4 FragColor;

uniform vec4 color;
uniform int isSurface;
uniform int useTexture;
uniform sampler2D diffuseTexture;

in VS_OUT {
  vec3 position;
  vec3 normal;
  vec2 texCoords;
  flat vec3 viewPosition;
} fs_in;

void main() {
  vec3 lightDirection = normalize(vec3(-1.0, 2.0, 0.0));
  vec3 normal = normalize(fs_in.normal);
  vec3 viewDirection = normalize(fs_in.viewPosition - fs_in.position);
  
  vec4 myColor = useTexture == 1 ? texture(diffuseTexture, fs_in.texCoords) : color;

  vec3 halfwayDirection = normalize(lightDirection + viewDirection);
  float ambient = 0.1;
  float normalDotLight = dot(normal, lightDirection);
  float diffuse = max(normalDotLight, 0.0);
  float specular = 0.2 * pow(max(dot(normal, halfwayDirection), 0.0), 32.0);
  float lighting = isSurface == 1 ? (ambient + diffuse + specular) : 1.0;
  FragColor = vec4(lighting * myColor.rgb, myColor.a);
}
