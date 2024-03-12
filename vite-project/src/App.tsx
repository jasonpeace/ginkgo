import "./App.css";
import React from "react";
import DNAQuery from "./components/DNAQuery";
import RequestDetail from "./components/RequestDetail";


import { BrowserRouter as Router, Route, Routes } from "react-router-dom";
//import { useCookies } from "react-cookie";

function App() {
  return (
    <Router>
      <Routes>
        <Route path="/detail/:id" element={<RequestDetail />} />
        <Route path="/" element={<DNAQuery />} />
      </Routes>
    </Router>
  );
}

export default App;
