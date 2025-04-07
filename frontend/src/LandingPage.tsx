import React from 'react';
import { motion } from 'framer-motion';
import { useNavigate } from 'react-router-dom';

const LandingPage: React.FC = () => {
  const navigate = useNavigate();

  return (
    <div className="min-h-screen bg-gradient-to-br from-green-50 to-blue-50 relative overflow-hidden">
      {/* Chemical structure background */}
      <div className="absolute inset-0 opacity-10">
        <div className="absolute top-0 left-0 w-full h-full">
          <div className="absolute top-1/4 left-1/4 w-32 h-32 border-4 border-green-300 rounded-full"></div>
          <div className="absolute top-1/3 right-1/4 w-24 h-24 border-4 border-blue-300 rounded-full"></div>
          <div className="absolute bottom-1/4 left-1/3 w-40 h-40 border-4 border-green-300 rounded-full"></div>
        </div>
      </div>

      {/* Main content */}
      <div className="relative z-10 max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-20">
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.8 }}
          className="text-center"
        >
          <h1 className="text-5xl md:text-6xl font-bold text-gray-900 mb-6">
            PolyMorph-AI
          </h1>
          <p className="text-xl md:text-2xl text-gray-700 mb-8">
            Revolutionizing Polymer Science with AI-Powered Sustainability
          </p>
          <p className="text-lg text-gray-600 max-w-3xl mx-auto mb-12">
            Our platform combines cutting-edge AI technology with polymer science to create sustainable materials for a better future. 
            Join us in transforming the way we design and manufacture polymers.
          </p>
          <motion.button
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
            onClick={() => navigate('/signup')}
            className="bg-green-600 hover:bg-green-700 text-white font-semibold py-3 px-8 rounded-lg text-lg transition-colors duration-300"
          >
            Get Started
          </motion.button>
        </motion.div>

        {/* Features section */}
        <div className="mt-24 grid grid-cols-1 md:grid-cols-3 gap-8">
          <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.8, delay: 0.2 }}
            className="bg-white p-6 rounded-xl shadow-lg"
          >
            <h3 className="text-xl font-semibold text-gray-900 mb-4">AI-Powered Analysis</h3>
            <p className="text-gray-600">Leverage advanced machine learning to predict polymer properties and optimize formulations.</p>
          </motion.div>
          <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.8, delay: 0.4 }}
            className="bg-white p-6 rounded-xl shadow-lg"
          >
            <h3 className="text-xl font-semibold text-gray-900 mb-4">Sustainable Solutions</h3>
            <p className="text-gray-600">Design eco-friendly polymers with reduced environmental impact and improved recyclability.</p>
          </motion.div>
          <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.8, delay: 0.6 }}
            className="bg-white p-6 rounded-xl shadow-lg"
          >
            <h3 className="text-xl font-semibold text-gray-900 mb-4">Collaborative Platform</h3>
            <p className="text-gray-600">Connect with researchers and industry experts to accelerate polymer innovation.</p>
          </motion.div>
        </div>
      </div>
    </div>
  );
};

export default LandingPage; 