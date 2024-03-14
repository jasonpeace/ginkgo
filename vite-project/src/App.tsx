import "./App.css";
import DNAQuery from "./components/DNAQuery";
import RequestDetail from "./components/RequestDetail";

import { BrowserRouter as Router, Route, Routes } from "react-router-dom";

function App() {

  const handleClick = () =>{
    location.href = "/";
  }

  return (
    <>
      <div className="header">
        <p className="header-title" onClick={handleClick}>Protein DNA Search</p>
      </div>
      <Router>
        <Routes>
          <Route path="/detail/:id" element={<RequestDetail />} />
          <Route path="/" element={<DNAQuery />} />
        </Routes>
      </Router>
    </>
  );
}

export default App;