import React from 'react';
import { useNavigate } from 'react-router-dom';

const LandingPage: React.FC = () => {
  const navigate = useNavigate();

  return (
    <div className="min-h-screen bg-gradient-to-br from-emerald-50 to-cyan-100 relative overflow-hidden">
      {/* faint grid background */}
      <div className="absolute inset-0 opacity-20 pointer-events-none">
        <svg className="w-full h-full" xmlns="http://www.w3.org/2000/svg">
          <defs>
            <pattern
              id="grid"
              x="0"
              y="0"
              width="40"
              height="40"
              patternUnits="userSpaceOnUse"
            >
              <path
                d="M 40 0 L 0 0 0 40"
                fill="none"
                stroke="#374151"
                strokeWidth="0.5"
              />
            </pattern>
          </defs>
          <rect width="100%" height="100%" fill="url(#grid)" />
        </svg>
      </div>

      {/* main content */}
      <div className="relative z-10 max-w-6xl mx-auto px-6 py-20 text-center">
        <h1 className="text-5xl md:text-6xl font-bold text-gray-900 mb-6">
          PolyMorph-AI
        </h1>
        <p className="text-xl md:text-2xl text-gray-700 mb-6">
          Self‑Decomposing Polymers Designed by Artificial Intelligence
        </p>
        <p className="text-lg text-gray-600 max-w-3xl mx-auto mb-10">
          Discover a new era of sustainable materials. Our AI‑driven platform
          creates polymers embedded with smart kill switches — activated by
          environmental cues like light, heat, or humidity — enabling safe,
          controlled decomposition and recyclability.
        </p>
        <button
          onClick={() => navigate('/signup')}
          className="bg-green-600 hover:bg-green-700 text-white font-semibold py-3 px-8 rounded-lg text-lg shadow-lg hover:shadow-xl transform hover:-translate-y-0.5 transition-all"
        >
          Get Started
        </button>

        {/* features */}
        <div className="mt-24 grid grid-cols-1 md:grid-cols-3 gap-8">
          <div className="bg-white bg-opacity-90 p-6 rounded-xl shadow-lg hover:shadow-xl transition-shadow">
            <h3 className="text-xl font-semibold text-gray-900 mb-4">
              AI‑Designed Molecules
            </h3>
            <p className="text-gray-600">
              Harness deep learning to generate polymer structures optimized for
              reactivity and eco‑safety.
            </p>
          </div>
          <div className="bg-white bg-opacity-90 p-6 rounded-xl shadow-lg hover:shadow-xl transition-shadow">
            <h3 className="text-xl font-semibold text-gray-900 mb-4">
              Kill Switch Engineering
            </h3>
            <p className="text-gray-600">
              Embed molecular triggers that activate under specific environmental
              conditions to initiate breakdown.
            </p>
          </div>
          <div className="bg-white bg-opacity-90 p-6 rounded-xl shadow-lg hover:shadow-xl transition-shadow">
            <h3 className="text-xl font-semibold text-gray-900 mb-4">
              Eco‑Conscious Design
            </h3>
            <p className="text-gray-600">
              Ensure lifecycle transparency and minimal waste through controlled
              material degradation.
            </p>
          </div>
        </div>
      </div>
    </div>
  );
};

export default LandingPage;
