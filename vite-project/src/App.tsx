import { useState, useEffect } from "react";
import "./App.css";
import "bootstrap/dist/css/bootstrap.min.css";
import { Button } from "react-bootstrap";
import Stack from "react-bootstrap/Stack";
import Table from "react-bootstrap/Table";

function App() {
  const [newDNAQuery, setNewDNAQuery] = useState();
  const [alignmentRequests, setAlignmentRequests] = useState([]);

  useEffect(() => {
    fetch_all_alignment_requests();
  }, []);

  function fetch_all_alignment_requests() {
    fetch("http://127.0.0.1:8000/requests")
      .then((response) => response.json())
      .then((data) => {
        console.log(data);
        setAlignmentRequests(data);
      });
  }

  function handleNewDNAQuery(event) {
    const data = event.target.value;
    setNewDNAQuery(data);
  }

  const handleSubmit = () => {
    fetch("http://127.0.0.1:8000/request", {
      method: "POST",
      headers: {
        Accept: "application/json",
        "Content-Type": "application/json",
      },
      body: JSON.stringify({ dna_string: newDNAQuery }),
    })
      .then(function (res) {
        fetch_all_alignment_requests();
        console.log(res);
      })
      .catch(function (res) {
        console.log(res);
      });
  };

  return (
    <>
      <Stack gap={3}>
        <label>Please enter your DNA query in the text field below.</label>
        <textarea
          name="DNAQueryForm"
          onChange={handleNewDNAQuery}
          value={newDNAQuery}
          className="dna-query"
        />
        <Button onClick={handleSubmit}>Submit DNA</Button>
      </Stack>

      {alignmentRequests.length > 0 && (
        <Table striped bordered hover>
          <thead>
            <tr>
              <td>ID</td>
              <td>DNA String</td>
              <td>Status</td>
              <td>Date Submitted</td>
              <td>Date Updated</td>
            </tr>
          </thead>
          <tbody>
            {alignmentRequests.map((request) => {
              return (
                <tr>
                  <td>{request.id}</td>
                  <td>{request.dna_string}</td>
                  <td>{request.status}</td>
                  <td>{request.date_submitted}</td>
                  <td>{request.date_updated}</td>
                </tr>
              );
            })}
          </tbody>
        </Table>
      )}
    </>
  );
}

export default App;
