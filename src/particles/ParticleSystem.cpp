//! \file ParticleSystem.cpp
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#include <GL/glew.h>
#include "../confparser/ConfParser.h"
#include "ParticleSystem.h"
#include "../lbm/LBM_def.h"


namespace particles {

ParticleSystem::ParticleSystem ( std::string configFileName )
    : numParticles_( 0 ),
      sizeBase_( 0. ),
      sizeVar_( 0. ),
      numSprites_( 0 ){
  try {

    ConfParser p;
    ConfBlock base = p.parse( configFileName );
    std::cout << "Parsed configuration file " << configFileName << std::endl;

    setup( base );
    solver_.setup( base );

    srand ( time(NULL) );

  } catch ( std::exception& e ) {
    std::cerr << e.what() << std::endl;
    exit( -1 );
  } catch ( const char* e ) {
    std::cerr << e << std::endl;
    exit( -1 );
  }

    minPart = vec3(-FLT_MAX,-FLT_MAX,-FLT_MAX);
    maxPart = vec3(FLT_MAX,FLT_MAX,FLT_MAX);
}

void ParticleSystem::setup( ConfBlock& base ) {
  try {

    // Read parameters from config file

    std::cout << "Setting up ParticleSystem..." << std::endl;

    // Read domain specification
    ConfBlock* paramBlock = base.find( "domain" );
    if ( paramBlock == NULL ) throw "No domain size given.";
    sizeX_ = paramBlock->getParam<int>( "sizeX" );
    sizeY_ = paramBlock->getParam<int>( "sizeY" );
    sizeZ_ = paramBlock->getParam<int>( "sizeZ" );
    std::cout << "Read domain specification:" << std::endl;
    std::cout << "sizeX : " << sizeX_ << std::endl;
    std::cout << "sizeY : " << sizeY_ << std::endl;
    std::cout << "sizeZ : " << sizeZ_ << std::endl;

    // Read parameter specification
    paramBlock = base.find( "parameters" );
    if ( paramBlock == NULL ) throw "No parameter specification given.";
    alpha_ = paramBlock->getParam<float>( "alpha" );
    beta_ = paramBlock->getParam<float>( "beta" );
    k_ = paramBlock->getParam<float>( "k" );
    float g_x = paramBlock->getParam<float>( "g_x" );
    float g_y = paramBlock->getParam<float>( "g_y" );
    float g_z = paramBlock->getParam<float>( "g_z" );
    gravity_ = vec3( g_x, g_y, g_z );
    smokeTemp_ = paramBlock->getParam<float>( "smokeTemp" );
    ambTemp_ = paramBlock->getParam<float>( "ambTemp" );
    maxSteps_ = paramBlock->getParam<int>( "maxSteps" );
    std::cout << "Read parameter specification:" << std::endl;
    std::cout << "alpha (conservation coefficient)   : " << alpha_ << std::endl;
    std::cout << "beta (transferability coefficient) : " << beta_ << std::endl;
    std::cout << "k (thermal expansion coefficient)  : " << k_ << std::endl;
    std::cout << "Gravity unit vector                : <" << g_x << "," << g_y << "," << g_z << ">" << std::endl;
    std::cout << "Temperature threshold for smoke (K): " << smokeTemp_ << std::endl;
    std::cout << "Ambient temperature (K)            : " << ambTemp_ << std::endl;
    std::cout << "Number of steps                    : " << maxSteps_ << std::endl;

    // Precompute Gauss function for thermal diffusion
    int maxElem = glm::sqrt( float(sizeX_ * sizeX_ + sizeY_ * sizeY_ + sizeZ_ * sizeZ_) );
    gaussTable_.reserve( maxElem );
    float a = 1. / glm::sqrt( 2. * M_PI );
    for ( int i = 0; i < maxElem; ++i ) {
      gaussTable_.push_back( a * exp(-  i * i / 2. ) );
    }

    // Set up particle system
    paramBlock = base.find( "particles" );
    if ( paramBlock == NULL ) throw "No particle emitters specified.";
    std::cout << "Set up the particle system..." << std::endl;

    // Optionally write particle update times to file
    paramBlock->getParam<std::string>( "updateTimeChart", updFileName_ );

    ConfBlock::childIterPair cip = paramBlock->findAll( "emitter" );

    float maxTemp = 0;
    for ( ConfBlock::childIter it = cip.first; it != cip.second; ++it ) {

      ConfBlock b = it->second;

      float xStart = b.getParam<float>( "xStart" );
      float xEnd   = b.getParam<float>( "xEnd" );
      float yStart = b.getParam<float>( "yStart" );
      float yEnd   = b.getParam<float>( "yEnd" );
      float zStart = b.getParam<float>( "zStart" );
      float zEnd   = b.getParam<float>( "zEnd" );
      float temp   = b.getParam<float>( "temp" );
      if ( temp > maxTemp ) maxTemp = temp;
      int fuel    = b.getParam<int>( "fuel" );
      int emitThreshold = b.getParam<int>( "emitThreshold" );
      float fuelConsumption = b.getParam<float>( "fuelConsumption" );
      float lifetimeCoeff = b.getParam<float>( "lifetimeCoeff" );
      std::cout << "Read emitter specification:" << std::endl;
      std::cout << "xStart : " << xStart << std::endl;
      std::cout << "xEnd   : " << xEnd << std::endl;
      std::cout << "yStart : " << yStart << std::endl;
      std::cout << "yEnd   : " << yEnd << std::endl;
      std::cout << "zStart : " << zStart << std::endl;
      std::cout << "zEnd   : " << zEnd << std::endl;
      std::cout << "temp : " << temp << std::endl;
      std::cout << "fuel : " << fuel << std::endl;
      std::cout << "emitThreshold : " << emitThreshold << std::endl;
      std::cout << "fuelConsumption : " << fuelConsumption << std::endl;
      std::cout << "lifetimeCoeff : " << lifetimeCoeff << std::endl;

      emitters_.push_back( Emitter( vec3( xStart, yStart, zStart ),
                                    vec3( xEnd - xStart,
                                                     yEnd - yStart,
                                                     zEnd - zStart ),
                                    temp,
                                    fuel,
                                    emitThreshold,
                                    fuelConsumption,
                                    lifetimeCoeff ) );

    }

  } catch ( std::exception& e ) {
    std::cerr << e.what() << std::endl;
    exit( -1 );
  } catch ( const char* e ) {
    std::cerr << e << std::endl;
    exit( -1 );
  }

  std::cout << "ParticleSystem setup finished!" << std::endl;
}

    void ParticleSystem::InitSDL()
    {
        if(SDL_Init(SDL_INIT_VIDEO) < 0)
        {
            std::cout << "Erreur lors de l'initialisation de la SDL : " << SDL_GetError() << std::endl;
            SDL_Quit();

        }


        // Version d'OpenGL

        SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);


        // Double Buffer

        SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
        SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);


        // Création de la fenêtre

        fenetre = SDL_CreateWindow("Test SDL 2.0", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 800, 600, SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL);

        if(fenetre == 0)
        {
            std::cout << "Erreur lors de la creation de la fenetre : " << SDL_GetError() << std::endl;
            SDL_Quit();
        }


        // Création du contexte OpenGL

        contexteOpenGL = SDL_GL_CreateContext(fenetre);

        if(contexteOpenGL == 0)
        {
            std::cout << SDL_GetError() << std::endl;
            SDL_DestroyWindow(fenetre);
            SDL_Quit();
        }

        GLenum initialisationGLEW( glewInit() );

        // Si l'initialisation a échouée :

        if(initialisationGLEW != GLEW_OK)
        {
            // On affiche l'erreur grâce à la fonction : glewGetErrorString(GLenum code)

            std::cout << "Erreur d'initialisation de GLEW : " << glewGetErrorString(initialisationGLEW) << std::endl;


            // On quitte la SDL

            SDL_GL_DeleteContext(contexteOpenGL);
            SDL_DestroyWindow(fenetre);
            SDL_Quit();
        }

        /* Enable smooth shading */
        glShadeModel( GL_SMOOTH );

        /* Set the background black */
        glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );

        /* Depth buffer setup */
        glClearDepth( 1.0f );

        /* Enables Depth Testing */
        glEnable( GL_DEPTH_TEST );

        /* The Type Of Depth Test To Do */
        glDepthFunc( GL_LEQUAL );

        /* Really Nice Perspective Calculations */
        glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );

        GLfloat ratio;

        int width = 800;
        int height = 600;

        /* Protect against a divide by zero */
        if ( height == 0 ) {
            height = 1;
        }

        ratio = ( GLfloat )width / ( GLfloat )height;

        /* Setup our viewport. */
        glViewport( 0, 0, ( GLsizei )width, ( GLsizei )height );

        /* change to the projection matrix and set our viewing volume. */
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity( );

        /* Set our perspective */
        gluPerspective( 45.0f, ratio, 0.1f, 100.0f );

        /* Make sure we're chaning the model view and not the projection */
        glMatrixMode( GL_MODELVIEW );

        /* Reset The View */
        glLoadIdentity( );
    }

    void ParticleSystem::CloseSDL()
    {
        SDL_GL_DeleteContext(contexteOpenGL);
        SDL_DestroyWindow(fenetre);
        SDL_Quit();
    }

    void ParticleSystem::Draw()
    {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if(event.type == SDL_QUIT)
                exit(0);
        }

        glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );
        /* Clear The Screen And The Depth Buffer */
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

        /* Move Left 1.5 Units And Into The Screen 6.0 */
        glLoadIdentity();

        gluLookAt(10,5,10,0,0,0,0,1,0);

        glPointSize(5);

        std::vector< Emitter >::iterator ite;
        std::list< Particle >::iterator itp;


        glBegin(GL_POINTS);
        glColor3f(1.0,1.0,1.0);
        for ( ite = emitters_.begin(); ite != emitters_.end(); ++ite ) {
            // Precompute modificator for particle size

            for (itp = (*ite).particles_.begin(); itp != (*ite).particles_.end(); ++itp) {
                Particle temp = (*itp);
                glVertex3f(temp.pos_.x,temp.pos_.y,temp.pos_.z);
                //std::cout << temp.pos_.x<< " " << temp.pos_.y << " " << temp.pos_.z << std::endl;
            }
        }
        glEnd();

        SDL_GL_SwapWindow(fenetre);

    }

void ParticleSystem::run() {

  // Discard first 20 steps to initialize velocity field
  for ( int i = 0; i < 20; ++i ) solver_.runStep();
  int step = 20;

    InitSDL();

  {

    float totTime = 0.;
    std::ofstream updFile;
    if ( updFileName_.length() ) {
      updFile.open( updFileName_.c_str(), std::ios::out );
      updFile << "\"Number of Particles\" \"Step time\" \"Emission time\" \"LBM step time\" \"Update time\"\n";
    }
    // Start the simulation loop
    for ( int step = 20; step < maxSteps_; ++step ) {

      struct timeval start, end;
      float stepTime = 0.;

      gettimeofday(&start, NULL);

      // emit particles
      emitParticles();

      gettimeofday(&end, NULL);
      float eTime = getTime( start, end );

      // Simulate one LBM step
      float sTime = solver_.runStep();

      gettimeofday(&start, NULL);

      // Update the particles
      updateParticles();

      gettimeofday(&end, NULL);
      float uTime = getTime( start, end );

      stepTime = eTime + sTime + uTime;

      std::cout << "Time step " << step << " of " << maxSteps_ << " took ";
      std::cout << eTime << " + " << sTime << " + " << uTime << " = ";
      std::cout << stepTime << "secs with " << numParticles_ << " particles";
      std::cout << std::endl;


        Draw();


      if ( updFileName_.length() )
        updFile << numParticles_ << " " << stepTime << " " << eTime << " " << sTime << " " << uTime << "\n";
    }

    if ( updFileName_.length() ) updFile.close();
  }

    CloseSDL();
    std::cout << "Particule min : " << minPart.x << " " << minPart.y << " " << minPart.z << std::endl;
    std::cout << "Particule max : " << maxPart.x << " " << maxPart.y << " " << maxPart.z << std::endl;
}

inline void ParticleSystem::updateParticles() {

  // Go over all particles of all emitters
  std::vector< Emitter >::iterator ite, ite2;
  std::list< Particle >::iterator itp, itp2;
  for ( ite = emitters_.begin(); ite != emitters_.end(); ++ite ) {
    // Precompute modificator for particle size
    float szCoeff = 1. / ( (*ite).fuel_ * (*ite).lifetimeCoeff_ );
    for ( itp = (*ite).particles_.begin(); itp != (*ite).particles_.end(); ++itp) {

        std::cout << "Begin ";
        std::cout << (*itp).pos_.x<< " " << (*itp).pos_.y << " " << (*itp).pos_.z << " Velo " ;

      // Remove particles that have left the domain or exceeded their lifetime
      while ( itp != (*ite).particles_.end() && ( (*itp).lifetime_ < 1 ||
              (*itp).getPos().x < 1 || (*itp).getPos().x > sizeX_ - 1 ||
              (*itp).getPos().y < 1 || (*itp).getPos().y > sizeY_ - 1 ||
              (*itp).getPos().z < 1 || (*itp).getPos().z > sizeZ_ -1  ||
              (*itp).temp_ < ambTemp_ ) ) {
//         std::cout << "Particle at pos <" << (*itp).getPos().X << "," << (*itp).getPos().Y << "," << (*itp).getPos().Z <<"> with lifetime " << (*itp).lifetime_ << " removed." << std::endl;

        itp = (*ite).particles_.erase( itp );
        numParticles_--;

          std::cout << " test ";
      }

      // If no particles are left anymore, go to next emitter
      if ( itp == (*ite).particles_.end() ) break;

      // Move particle according to fluid velocity at current position
      Vec3<float> vel = solver_.getVelocity( (*itp).pos_.x, (*itp).pos_.y, (*itp).pos_.z );
      vec3 v = vec3( vel[0], vel[1], vel[2] );


      // Buoyancy force
      vec3 g = -gravity_ * k_ * ( (*itp).temp_ - ambTemp_ );
      // Check whether the buoyancy force would carry the particle into an
      // obstacle

        std::cout << g.x<< " " << g.y << " " << g.z << " Velo " ;

      (*itp).updatePos( v + g );

      // Update temperature
      float tempExt = - gaussTable_[0] * (*itp).temp_;
      // Add up temperature contributions of all other particles weighted by distance
      for ( ite2 = emitters_.begin(); ite2 != emitters_.end(); ++ite2 ) {
        for ( itp2 = (*ite2).particles_.begin(); itp2 != (*ite2).particles_.end(); ++itp2) {
          tempExt += gaussTable_[ (int) (*itp).dist( *itp2 ) ] * (*itp2).temp_;
        }
      }
      (*itp).temp_ = alpha_ * (*itp).temp_ + beta_ * tempExt;

      if ( (*itp).type_ == FIRE ) {

        // If temperature has fallen below threshold, convert to smoke particle
        if ( (*itp).temp_ < smokeTemp_ ) {
          (*itp).setSmoke( szCoeff );
        }

      }

      (*itp).lifetime_--;

        minPart.x = glm::min(minPart.x,(*itp).pos_.x);
        minPart.y = glm::min(minPart.y,(*itp).pos_.y);
        minPart.z = glm::min(minPart.z,(*itp).pos_.z);

        maxPart.x = glm::max(maxPart.x,(*itp).pos_.x);
        maxPart.y = glm::max(maxPart.y,(*itp).pos_.y);
        maxPart.z = glm::max(maxPart.z,(*itp).pos_.z);

        std::cout << (*itp).pos_.x<< " " << (*itp).pos_.y << " " << (*itp).pos_.z << " end " << std::endl;
    }
  }

  // Move the spheres
  for ( uint i = 0; i < spheres_.size(); ++i ) {
    spheres_[i].move();
  }

}

inline void ParticleSystem::emitParticles() {
  // Go over all emitters
  std::vector< Emitter >::iterator ite;
  for ( ite = emitters_.begin(); ite != emitters_.end(); ++ite ) {
    // Emit random number of particles up to emit threshold
    for ( int i = 0; i < (*ite).emitThreshold_; ++i ) {
      vec3 pos = (*ite).pos_ + vec3(
          ( std::rand() * (*ite).size_.x ) / (float) RAND_MAX,
          ( std::rand() * (*ite).size_.y ) / (float) RAND_MAX,
          ( std::rand() * (*ite).size_.z ) / (float) RAND_MAX
                                                         );



        (*ite).particles_.push_back(
            Particle( pos, (*ite).temp_, (int) ((*ite).fuel_ * (*ite).lifetimeCoeff_ ) ) );
      // Reduce emitter's fuel
      (*ite).fuel_ *= (*ite).fuelConsumption_;
      numParticles_++;
    }
  }
}

void ParticleSystem::generateBlackBodyColorTable( float maxTemp ) {

  float cie_colour_match[81][3] = {
      {0.0014,0.0000,0.0065}, {0.0022,0.0001,0.0105}, {0.0042,0.0001,0.0201},
      {0.0076,0.0002,0.0362}, {0.0143,0.0004,0.0679}, {0.0232,0.0006,0.1102},
      {0.0435,0.0012,0.2074}, {0.0776,0.0022,0.3713}, {0.1344,0.0040,0.6456},
      {0.2148,0.0073,1.0391}, {0.2839,0.0116,1.3856}, {0.3285,0.0168,1.6230},
      {0.3483,0.0230,1.7471}, {0.3481,0.0298,1.7826}, {0.3362,0.0380,1.7721},
      {0.3187,0.0480,1.7441}, {0.2908,0.0600,1.6692}, {0.2511,0.0739,1.5281},
      {0.1954,0.0910,1.2876}, {0.1421,0.1126,1.0419}, {0.0956,0.1390,0.8130},
      {0.0580,0.1693,0.6162}, {0.0320,0.2080,0.4652}, {0.0147,0.2586,0.3533},
      {0.0049,0.3230,0.2720}, {0.0024,0.4073,0.2123}, {0.0093,0.5030,0.1582},
      {0.0291,0.6082,0.1117}, {0.0633,0.7100,0.0782}, {0.1096,0.7932,0.0573},
      {0.1655,0.8620,0.0422}, {0.2257,0.9149,0.0298}, {0.2904,0.9540,0.0203},
      {0.3597,0.9803,0.0134}, {0.4334,0.9950,0.0087}, {0.5121,1.0000,0.0057},
      {0.5945,0.9950,0.0039}, {0.6784,0.9786,0.0027}, {0.7621,0.9520,0.0021},
      {0.8425,0.9154,0.0018}, {0.9163,0.8700,0.0017}, {0.9786,0.8163,0.0014},
      {1.0263,0.7570,0.0011}, {1.0567,0.6949,0.0010}, {1.0622,0.6310,0.0008},
      {1.0456,0.5668,0.0006}, {1.0026,0.5030,0.0003}, {0.9384,0.4412,0.0002},
      {0.8544,0.3810,0.0002}, {0.7514,0.3210,0.0001}, {0.6424,0.2650,0.0000},
      {0.5419,0.2170,0.0000}, {0.4479,0.1750,0.0000}, {0.3608,0.1382,0.0000},
      {0.2835,0.1070,0.0000}, {0.2187,0.0816,0.0000}, {0.1649,0.0610,0.0000},
      {0.1212,0.0446,0.0000}, {0.0874,0.0320,0.0000}, {0.0636,0.0232,0.0000},
      {0.0468,0.0170,0.0000}, {0.0329,0.0119,0.0000}, {0.0227,0.0082,0.0000},
      {0.0158,0.0057,0.0000}, {0.0114,0.0041,0.0000}, {0.0081,0.0029,0.0000},
      {0.0058,0.0021,0.0000}, {0.0041,0.0015,0.0000}, {0.0029,0.0010,0.0000},
      {0.0020,0.0007,0.0000}, {0.0014,0.0005,0.0000}, {0.0010,0.0004,0.0000},
      {0.0007,0.0002,0.0000}, {0.0005,0.0002,0.0000}, {0.0003,0.0001,0.0000},
      {0.0002,0.0001,0.0000}, {0.0002,0.0001,0.0000}, {0.0001,0.0000,0.0000},
      {0.0001,0.0000,0.0000}, {0.0001,0.0000,0.0000}, {0.0000,0.0000,0.0000}
  };

  for ( float t = smokeTemp_; t < maxTemp; t += 50. ) {

    // Calculate x,y,z colors from solar spectrum

    int i;
    float lambda, x = 0, y = 0, z = 0, xyz;
    for (i = 0, lambda = 380; lambda < 780.1; i++, lambda += 5) {
        float Me;
        // Get black body radiation intensity for given temperature and
        // wavelength
        float wlm = lambda * 1e-9; // wavelength in meters
        Me = (3.74183e-16 * pow(wlm, -5.0)) / (exp(1.4388e-2 / (wlm * t)) - 1.0);
        x += Me * cie_colour_match[i][0];
        y += Me * cie_colour_match[i][1];
        z += Me * cie_colour_match[i][2];
    }
    xyz = (x + y + z);
    x /= xyz;
    y /= xyz;
    z /= xyz;

    // Calculate r,g,b colors from x,y,z colors

    float xr, yr, zr, xg, yg, zg, xb, yb, zb;
    float xw, yw, zw;
    float rx, ry, rz, gx, gy, gz, bx, by, bz;
    float rw, gw, bw;

    xr = 0.630;  yr = 0.340;  zr = 1 - (xr + yr);
    xg = 0.310;  yg = 0.595;  zg = 1 - (xg + yg);
    xb = 0.155;  yb = 0.070;  zb = 1 - (xb + yb);

    xw = 0.3127; yw = 0.3291; zw = 1 - (xw + yw);

    // xyz -> rgb matrix, before scaling to white.
    rx = (yg * zb) - (yb * zg);  ry = (xb * zg) - (xg * zb);  rz = (xg * yb) - (xb * yg);
    gx = (yb * zr) - (yr * zb);  gy = (xr * zb) - (xb * zr);  gz = (xb * yr) - (xr * yb);
    bx = (yr * zg) - (yg * zr);  by = (xg * zr) - (xr * zg);  bz = (xr * yg) - (xg * yr);

    // White scaling factors.
    // Dividing by yw scales the white luminance to unity, as conventional.
    rw = ((rx * xw) + (ry * yw) + (rz * zw)) / yw;
    gw = ((gx * xw) + (gy * yw) + (gz * zw)) / yw;
    bw = ((bx * xw) + (by * yw) + (bz * zw)) / yw;

    // xyz -> rgb matrix, correctly scaled to white.
    rx = rx / rw;  ry = ry / rw;  rz = rz / rw;
    gx = gx / gw;  gy = gy / gw;  gz = gz / gw;
    bx = bx / bw;  by = by / bw;  bz = bz / bw;

    // rgb of the desired point
    float r = (rx * x) + (ry * y) + (rz * z);
    float g = (gx * x) + (gy * y) + (gz * z);
    float b = (bx * x) + (by * y) + (bz * z);

    // Constrain to RGB gammut
    double w;
    // Amount of white needed is w = - min(0, *r, *g, *b)
    w = (0 < r) ? 0 : r;
    w = (w < g) ? w : g;
    w = (w < b) ? w : b;
    w = -w;
    // Add just enough white to make r, g, b all positive.
    if (w > 0) {
        r += w;  g += w; b += w;
    }

    // Normalize
    float max = (r > b) ? ( (r > g) ? r : g ) : ( (b > g) ? b : g );
    r *= 255. / max;
    g *= 255. / max;
    b *= 255. / max;

    std::cout << "Temperature " << t << "K: RGB <" << r << ", " << g << ", " << b << ">" << std::endl;
  }
  std::cout << "Generated black body color table from " << smokeTemp_;
  std::cout << "K to " << maxTemp << "K (" <<  " values)" << std::endl;
}

} // namespace particles
