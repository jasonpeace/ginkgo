import { useState, useEffect } from "react";
import { useParams } from 'react-router-dom';
import Table from "react-bootstrap/Table";
import React from "react";

function RequestDetail() {
  let { id } = useParams();
  const [details, setDetails] = useState();
  const [result, setResult] = useState();
  const [detailId, setDetailId] = useState(id);
  
  useEffect(() => {
    fetch_details();
  }, [detailId]);

  function fetch_details() {
    fetch("http://127.0.0.1:8000/api/detail/" + detailId, {
      method: "GET",
    })
      .then((response) => response.json())
      .then((data) => {
        console.log(data);
        setDetails(data.alignment_request);
      });
  }

  function handleUpdateRequest(data) {
    setRequest(data);
    console.log(data);
  }

  function handleUpdateResult(data) {
    setResult(data);
    console.log(data);
  }

  return (
    <>
      {details && (
        <Table striped bordered hover>
          <thead>
            <tr>
              <td>DNA String</td>
              <td>Status</td>
              <td>Date Submitted</td>
              <td>Date Updated</td>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td>{details.dna_string}</td>
              <td>{details.status}</td>
              <td>{details.date_submitted}</td>
              <td>{details.date_updated}</td>
            </tr>
          </tbody>
        </Table>
      )}

      {/* <Table striped bordered hover>
        <thead>
          <tr>
            <td>Protein id</td>
            <td>Alignment Result</td>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td>{details.protein_id}</td>
            <td>{details.alignment_result}</td>
          </tr>
        </tbody>
      </Table> */}
    </>
  );
}

export default RequestDetail;
