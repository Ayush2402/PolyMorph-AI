import React, { useState, useEffect, useRef } from 'react';
import axios from 'axios';
import { motion, AnimatePresence } from 'framer-motion';

declare global {
  interface Window {
    SmilesDrawer: any;
  }
}

interface PolymerResponse {
  smiles: string;
  has_kill_switch: boolean;
  degradability_score: number;
  nomenclature?: {
    name: string;
    formula: string;
    functional_groups: string[];
    structure_type: string;
  };
  visualization?: {
    atoms: Array<{
      symbol: string;
      atomic_num: number;
      formal_charge: number;
      is_aromatic: boolean;
      hybridization: string;
    }>;
    bonds: Array<{
      bond_type: string;
      is_conjugated: boolean;
      is_aromatic: boolean;
      begin_atom_idx: number;
      end_atom_idx: number;
    }>;
    num_rings: number;
    molecular_weight: number;
    num_rotatable_bonds: number;
    num_h_donors: number;
    num_h_acceptors: number;
  };
}

interface TriggerOption {
  value: string;
  label: string;
}

function MainApp() {
  // State variables
  const [domain, setDomain] = useState('agriculture');
  const [trigger, setTrigger] = useState('');
  const [result, setResult] = useState<PolymerResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [activeTab, setActiveTab] = useState('atoms'); // 'atoms' or 'bonds'
  const canvasRef = useRef<HTMLCanvasElement>(null);

  // Domain-specific trigger options
  const triggerOptions: Record<string, TriggerOption[]> = {
    agriculture: [
      { value: 'pH < 5', label: 'Low pH (< 5)' },
      { value: 'temperature > 30', label: 'High Temperature (> 30°C)' },
      { value: 'microbial_presence', label: 'Microbial Presence' }
    ],
    water: [
      { value: 'pH > 8', label: 'High pH (> 8)' },
      { value: 'UV_light', label: 'UV Light Exposure' },
      { value: 'oxidation', label: 'Oxidative Environment' }
    ],
    urban: [
      { value: 'temperature > 40', label: 'High Temperature (> 40°C)' },
      { value: 'mechanical_stress', label: 'Mechanical Stress' },
      { value: 'moisture', label: 'Moisture Exposure' }
    ]
  };

  const handleGenerate = async () => {
    try {
      setLoading(true);
      setError('');
      const response = await axios.post('http://localhost:8000/api/generate', {
        domain,
        trigger_condition: trigger || undefined
      });
      setResult(response.data);
    } catch (err: any) {
      console.error('Error details:', err);
      if (err.response) {
        // The request was made and the server responded with a status code
        // that falls out of the range of 2xx
        setError(`Server error: ${err.response.data.detail || 'Unknown error'}`);
      } else if (err.request) {
        // The request was made but no response was received
        setError('No response from server. Please make sure the backend is running.');
      } else {
        // Something happened in setting up the request that triggered an Error
        setError(`Error: ${err.message}`);
      }
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    if (result?.smiles && canvasRef.current) {
      try {
        const SmilesDrawer = window.SmilesDrawer;
        if (!SmilesDrawer) {
          console.error('SmilesDrawer library not loaded');
          return;
        }

        // Initialize the drawer
        const drawer = new SmilesDrawer.Drawer({ width: 500, height: 300 });
        // Parse the SMILES string
        SmilesDrawer.parse(result.smiles, function(tree: any) {
          // Draw the parsed structure
          drawer.draw(tree, canvasRef.current, 'light', false);
        });
      } catch (err) {
        console.error('Error drawing molecule:', err);
      }
    }
  }, [result?.smiles]);

  const renderMoleculeInfo = () => {
    if (!result?.visualization) return null;

    const { atoms, bonds, num_rings, molecular_weight, num_rotatable_bonds, num_h_donors, num_h_acceptors } = result.visualization;
    
    return (
      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.5 }}
        className="p-6 bg-white rounded-xl shadow-md"
      >
        <motion.h3 
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={{ delay: 0.2 }}
          className="text-xl font-semibold mb-4"
        >
          Molecule Information
        </motion.h3>
        
        {/* Nomenclature Section */}
        {result.nomenclature && (
          <motion.div 
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ delay: 0.4, duration: 0.3 }}
            className="mb-6 p-4 bg-green-50 rounded-lg"
          >
            <h4 className="font-medium mb-3">Nomenclature</h4>
            <div className="grid grid-cols-1 gap-2 text-sm">
              <p><span className="font-medium">Name:</span> {result.nomenclature.name}</p>
              <p><span className="font-medium">Formula:</span> {result.nomenclature.formula}</p>
              <p><span className="font-medium">Structure Type:</span> {result.nomenclature.structure_type}</p>
              <p>
                <span className="font-medium">Functional Groups:</span>{' '}
                {result.nomenclature.functional_groups.length > 0 
                  ? result.nomenclature.functional_groups.join(', ')
                  : 'None identified'}
              </p>
            </div>
          </motion.div>
        )}
        
        {/* Structure Visualization */}
        <motion.div 
          initial={{ opacity: 0, scale: 0.95 }}
          animate={{ opacity: 1, scale: 1 }}
          transition={{ delay: 0.3, duration: 0.5 }}
          className="mb-6"
        >
          <h4 className="font-medium mb-2">Structural Visualization</h4>
          <canvas 
            ref={canvasRef} 
            className="border border-gray-300 mx-auto bg-white rounded-lg shadow-sm"
            width="500"
            height="300"
          ></canvas>
        </motion.div>
        
        {/* Molecular Properties */}
        <motion.div 
          initial={{ opacity: 0, y: 10 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ delay: 0.5, duration: 0.3 }}
          className="mb-6 p-4 bg-gray-50 rounded-lg"
        >
          <h4 className="font-medium mb-3">Molecular Properties</h4>
          <div className="grid grid-cols-2 gap-4">
            <div className="text-sm">
              <p><span className="font-medium">Molecular Weight:</span> {molecular_weight.toFixed(2)} g/mol</p>
              <p><span className="font-medium">Number of Rings:</span> {num_rings}</p>
              <p><span className="font-medium">Rotatable Bonds:</span> {num_rotatable_bonds}</p>
            </div>
            <div className="text-sm">
              <p><span className="font-medium">H-bond Donors:</span> {num_h_donors}</p>
              <p><span className="font-medium">H-bond Acceptors:</span> {num_h_acceptors}</p>
            </div>
          </div>
        </motion.div>

        {/* Tabs for Atoms and Bonds */}
        <div className="mb-4">
          <div className="border-b border-gray-200">
            <nav className="flex -mb-px">
              <button
                onClick={() => setActiveTab('atoms')}
                className={`py-2 px-4 text-center border-b-2 font-medium text-sm ${
                  activeTab === 'atoms'
                    ? 'border-blue-500 text-blue-600'
                    : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
                }`}
              >
                Atoms ({atoms.length})
              </button>
              <button
                onClick={() => setActiveTab('bonds')}
                className={`py-2 px-4 text-center border-b-2 font-medium text-sm ${
                  activeTab === 'bonds'
                    ? 'border-blue-500 text-blue-600'
                    : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
                }`}
              >
                Bonds ({bonds.length})
              </button>
            </nav>
          </div>
          
          <AnimatePresence mode="wait">
            {activeTab === 'atoms' ? (
              <motion.div
                key="atoms"
                initial={{ opacity: 0, y: 10 }}
                animate={{ opacity: 1, y: 0 }}
                exit={{ opacity: 0, y: -10 }}
                transition={{ duration: 0.2 }}
                className="mt-4"
              >
                <div className="space-y-1 max-h-60 overflow-y-auto p-2">
                  {atoms.map((atom, idx) => (
                    <motion.div 
                      key={idx} 
                      initial={{ opacity: 0, x: -5 }}
                      animate={{ opacity: 1, x: 0 }}
                      transition={{ delay: idx * 0.03, duration: 0.2 }}
                      className="text-sm p-2 border-b border-gray-100 last:border-b-0"
                    >
                      {atom.symbol} ({atom.atomic_num}) - {atom.hybridization}
                      {atom.is_aromatic && ' (Aromatic)'}
                      {atom.formal_charge !== 0 && ` [${atom.formal_charge}]`}
                    </motion.div>
                  ))}
                </div>
              </motion.div>
            ) : (
              <motion.div
                key="bonds"
                initial={{ opacity: 0, y: 10 }}
                animate={{ opacity: 1, y: 0 }}
                exit={{ opacity: 0, y: -10 }}
                transition={{ duration: 0.2 }}
                className="mt-4"
              >
                <div className="space-y-1 max-h-60 overflow-y-auto p-2">
                  {bonds.map((bond, idx) => (
                    <motion.div 
                      key={idx} 
                      initial={{ opacity: 0, x: -5 }}
                      animate={{ opacity: 1, x: 0 }}
                      transition={{ delay: idx * 0.03, duration: 0.2 }}
                      className="text-sm p-2 border-b border-gray-100 last:border-b-0"
                    >
                      {bond.bond_type} between atoms {bond.begin_atom_idx} and {bond.end_atom_idx}
                      {bond.is_aromatic && ' (Aromatic)'}
                      {bond.is_conjugated && ' (Conjugated)'}
                    </motion.div>
                  ))}
                </div>
              </motion.div>
            )}
          </AnimatePresence>
        </div>
      </motion.div>
    );
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100">
      <div className="container mx-auto px-4 py-8">
        <motion.h1
          initial={{ opacity: 0, y: -20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5 }}
          className="text-3xl font-bold text-center text-gray-800 mb-8"
        >
          AI Self‑Decomposing Polymers
        </motion.h1>
        
        <div className="flex flex-col lg:flex-row gap-6">
          {/* Input Panel */}
          <motion.div
            initial={{ opacity: 0, x: -20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ duration: 0.5, delay: 0.2 }}
            className="lg:w-1/2 bg-white p-6 rounded-xl shadow-md"
          >
            <h2 className="text-xl font-semibold mb-6">Polymer Generation</h2>
            
            <div className="mb-6">
              <label className="block text-gray-700 text-sm font-bold mb-2">
                Domain
              </label>
              <select
                value={domain}
                onChange={(e) => {
                  setDomain(e.target.value);
                  setTrigger(''); // Reset trigger when domain changes
                }}
                className="shadow border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-transparent"
              >
                <option value="agriculture">Agriculture</option>
                <option value="water">Water</option>
                <option value="urban">Urban</option>
              </select>
            </div>

            <div className="mb-6">
              <label className="block text-gray-700 text-sm font-bold mb-2">
                Trigger Condition
              </label>
              <div className="space-y-2">
                <div className="flex flex-wrap gap-2">
                  {triggerOptions[domain].map((option) => (
                    <motion.button
                      key={option.value}
                      type="button"
                      whileHover={{ scale: 1.05 }}
                      whileTap={{ scale: 0.95 }}
                      onClick={() => setTrigger(option.value)}
                      className={`px-3 py-2 text-sm rounded-lg transition-colors ${
                        trigger === option.value
                          ? 'bg-blue-500 text-white'
                          : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
                      }`}
                    >
                      {option.label}
                    </motion.button>
                  ))}
                  {trigger && (
                    <motion.button
                      initial={{ opacity: 0, scale: 0.9 }}
                      animate={{ opacity: 1, scale: 1 }}
                      whileHover={{ scale: 1.05 }}
                      whileTap={{ scale: 0.95 }}
                      onClick={() => setTrigger('')}
                      className="px-3 py-2 text-sm rounded-lg bg-red-100 text-red-700 hover:bg-red-200"
                    >
                      Clear
                    </motion.button>
                  )}
                </div>
                
                <div className="text-xs text-gray-500 italic">
                  Or specify a custom condition:
                </div>
                <input
                  type="text"
                  value={trigger}
                  onChange={(e) => setTrigger(e.target.value)}
                  className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                  placeholder="e.g., pH < 5, temperature > 30"
                />
              </div>
            </div>

            <motion.button
              whileHover={{ scale: 1.03 }}
              whileTap={{ scale: 0.97 }}
              disabled={loading}
              onClick={handleGenerate}
              className="w-full bg-gradient-to-r from-blue-500 to-blue-600 text-white font-bold py-3 px-4 rounded-lg focus:outline-none focus:shadow-outline transition-all duration-300 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              {loading ? (
                <div className="flex items-center justify-center">
                  <div className="w-5 h-5 border-2 border-white border-t-transparent rounded-full animate-spin mr-2"></div>
                  Generating...
                </div>
              ) : (
                'Generate Polymer'
              )}
            </motion.button>

            {error && (
              <motion.div
                initial={{ opacity: 0, y: 10 }}
                animate={{ opacity: 1, y: 0 }}
                className="text-red-500 text-sm mt-4 p-3 bg-red-50 rounded-lg"
              >
                {error}
              </motion.div>
            )}
          </motion.div>
          
          {/* Results Panel */}
          <motion.div
            initial={{ opacity: 0, x: 20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ duration: 0.5, delay: 0.3 }}
            className="lg:w-1/2"
          >
            {result ? (
              <motion.div
                layout
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                className="bg-white p-6 rounded-xl shadow-md mb-6"
              >
                <motion.h2 
                  initial={{ opacity: 0 }}
                  animate={{ opacity: 1 }}
                  transition={{ delay: 0.1 }}
                  className="text-xl font-semibold mb-4"
                >
                  Generated Polymer
                </motion.h2>
                <div className="space-y-4">
                  <motion.div
                    initial={{ opacity: 0, y: 10 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ delay: 0.2 }}
                    className="p-3 bg-gray-50 rounded-lg"
                  >
                    <div className="font-medium text-gray-700">SMILES</div>
                    <div className="font-mono text-sm break-all mt-1">{result.smiles}</div>
                  </motion.div>
                  
                  <div className="grid grid-cols-2 gap-4">
                    <motion.div
                      initial={{ opacity: 0, y: 10 }}
                      animate={{ opacity: 1, y: 0 }}
                      transition={{ delay: 0.3 }}
                      className="p-3 bg-gray-50 rounded-lg"
                    >
                      <div className="font-medium text-gray-700">Kill Switch</div>
                      <div className={`font-semibold mt-1 ${result.has_kill_switch ? 'text-green-600' : 'text-red-600'}`}>
                        {result.has_kill_switch ? 'Yes' : 'No'}
                      </div>
                    </motion.div>
                    
                    <motion.div
                      initial={{ opacity: 0, y: 10 }}
                      animate={{ opacity: 1, y: 0 }}
                      transition={{ delay: 0.4 }}
                      className="p-3 bg-gray-50 rounded-lg"
                    >
                      <div className="font-medium text-gray-700">Degradability Score</div>
                      <div className="font-semibold mt-1">
                        <span
                          className={
                            result.degradability_score > 0.7
                              ? 'text-green-600'
                              : result.degradability_score > 0.4
                              ? 'text-yellow-600'
                              : 'text-red-600'
                          }
                        >
                          {(result.degradability_score * 100).toFixed(1)}%
                        </span>
                      </div>
                    </motion.div>
                  </div>
                </div>
              </motion.div>
            ) : (
              <div className="bg-gray-50 p-6 rounded-xl border-2 border-dashed border-gray-300 flex flex-col items-center justify-center min-h-[300px]">
                <svg className="w-16 h-16 text-gray-400 mb-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
                </svg>
                <p className="text-gray-500 mb-2">No polymer generated yet</p>
                <p className="text-gray-400 text-sm text-center">
                  Configure your settings and click Generate Polymer
                </p>
              </div>
            )}
            
            {result?.visualization && renderMoleculeInfo()}
          </motion.div>
        </div>
      </div>
    </div>
  );
}

export default MainApp; 